#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import numpy as np
import pysam
import time

# Import VaSe specific classes.
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile
from VcfVariant import VcfVariant


class VaSeBuilder:
    # Constructor that saves the identifier, date and time of the current
    def __init__(self, vaseid):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.creation_id = str(vaseid)
        self.creation_date = datetime.now().date()
        self.creation_time = datetime.now().time()

        # Create the bookkeeping variables used for saving the variant contexts.
        # VariantContextFile that saves the acceptor, donor and variant contexts and their associated data.
        self.contexts = VariantContextFile()
        self.vaselogger.info(f"VaSeBuilder: {self.creation_id} ; {self.creation_date} ; {self.creation_time}")

        # Dictionary used for debug messages.
        self.debug_dict = {"vw": "Establishing search window",
                           "dr": "Gathering donor variant reads",
                           "dc": "Determining donor context",
                           "ar": "Gathering acceptor variant reads",
                           "ac": "Determining acceptor context",
                           "cc": "Determining combined context",
                           "cdr": "Gathering combined context donor reads",
                           "car": "Gathering combined context acceptor reads",
                           "done": "Variant complete. Processing"
                           }

    # Method to print debug messages.
    def debug_msg(self, step, variant_id, t0=False):
        process = self.debug_dict[step]
        for_var = f"for variant {variant_id}"
        if t0:
            took = f" took {time.time() - t0} seconds."
        else:
            took = '.'
        self.vaselogger.debug(f"{process} {for_var}{took}")

    # ===METHODS TO PRODUCE VARIANT CONTEXTS===================================
    # Creates the new FastQ validation dataset by replacing acceptor reads
    # containing a VCF variant with BAM reads from patients.  Returns
    # true at the end to indicate the process is done.

    def build_varcon_set(self, sampleidlist,
                         vcfsamplemap, bamsamplemap,
                         acceptorbamloc,
                         outpath,
                         reference_loc,
                         varcon_outpath,
                         variant_list):
        self.vaselogger.info("Begin building the validation set.")
        start_time = self.creation_time
        donor_vcfs_used, donor_bams_used = [], []

        try:
            acceptorbamfile = pysam.AlignmentFile(acceptorbamloc, reference_filename=reference_loc)
        except IOError:
            self.vaselogger.critical("Could not open acceptor BAM file. Exitting.")
            exit()

        # Iterate over the samples to use for building the validation set.
        for sampleid in sampleidlist:
            self.vaselogger.debug(f"Processing data for sample {sampleid}.")

            sample_variant_filter = None
            if variant_list is not None:
                if sampleid in variant_list:
                    sample_variant_filter = variant_list[sampleid]
            samplevariants = self.get_sample_vcf_variants(vcfsamplemap[sampleid], sample_variant_filter)
            if len(samplevariants) == 0:
                self.vaselogger.warning(f"No variants obtained for sample {sampleid}. Skipping sample")
                continue

            # Check in the VCF BAM link map that there is a BAM file for the sample as well.
            try:
                bamfile = pysam.AlignmentFile(bamsamplemap[sampleid], reference_filename=reference_loc)
                self.vaselogger.debug("Opened BAM file "
                                      f"{bamsamplemap[sampleid]}")
            # Skip this sample if there is a problem opening the BAM or VCF.
            except IOError:
                self.vaselogger.warning(f"Could not open data files for sample {sampleid}. Skipping sample.")
                continue

            # Loop over the variants in the VCF file.
            for vcfvar in samplevariants:
                var_t0 = time.time()

                variantid = vcfvar.get_variant_id()
                acceptor_unmapped = []
                donor_unmapped = []
                varcon_unmapped_a = []
                varcon_unmapped_d = []
                self.debug_msg("vw", variantid)

                # Determine the search window before gathering reads overlapping with the VCF variant.
                self.vaselogger.debug(f"Variant {vcfvar.get_variant_id()} determined to be a "
                                      f"{vcfvar.get_variant_type()}.")
                searchwindow = self.determine_read_search_window(vcfvar.get_variant_type(),
                                                                 vcfvar)
                self.vaselogger.debug(f"Search window determined to be {vcfvar.get_variant_chrom()}:{searchwindow[0]+1}"
                                      f"-{searchwindow[1]}.")

                # Check if variant is already in an in-use context;
                # If so, skip it.
                if (self.contexts.variant_is_in_context(vcfvar.get_variant_type(),
                                                        vcfvar.get_variant_chrom(),
                                                        searchwindow[0],
                                                        searchwindow[1])):
                    self.vaselogger.debug(f"VCF variant {variantid} is located in an already existing variant context.")
                    continue

                try:
                    # === DONOR ===============================================
                    # Obtain all donor BAM reads at the variant position, as well as their mates.
                    self.debug_msg("dr", variantid)
                    t0 = time.time()
                    donor_context_reads = (
                            self.get_variant_reads(variantid,
                                                   vcfvar.get_variant_chrom(),
                                                   searchwindow[0],
                                                   searchwindow[1],
                                                   bamfile,
                                                   True,
                                                   donor_unmapped)
                            )
                    self.debug_msg("dr", variantid, t0)

                    # If no donor reads were found at this position, then skip this variant.
                    if not donor_context_reads:
                        self.vaselogger.info(f"No reads found for variant {variantid} in donor {sampleid}. "
                                             f"Skipping variant.")
                        continue

                    # Determine the donor context based on the reads overlapping the variant.
                    self.debug_msg("dc", variantid)
                    t0 = time.time()
                    donor_context = (
                            self.determine_context(donor_context_reads,
                                                   vcfvar.get_variant_pos(),
                                                   vcfvar.get_variant_chrom())
                            )
                    self.debug_msg("dc", variantid, t0)

                    # === ACCEPTOR ============================================
                    # Obtain all acceptor BAM reads at the variant position, as well as their mates.
                    self.debug_msg("ar", variantid)
                    t0 = time.time()
                    acceptor_context_reads = (
                            self.get_variant_reads(variantid,
                                                   vcfvar.get_variant_chrom(),
                                                   searchwindow[0],
                                                   searchwindow[1],
                                                   acceptorbamfile,
                                                   True,
                                                   acceptor_unmapped)
                            )
                    self.debug_msg("ar", variantid, t0)

                    # Only perform the following if no reads were found in
                    # the acceptor at the variant location.
                    if not acceptor_context_reads:
                        self.vaselogger.warning(
                                f"No reads found for variant {variantid} in "
                                "acceptor. Acceptor and donor sequencing may "
                                "have been performed with different methods. "
                                "Proceeding anyway."
                                )
                        # Create a dummy list of reads.
                        acceptor_context_reads = None
                        # Temporarily set acceptor context equal to donor.
                        acceptor_context = donor_context

                    # If reads WERE found, proceed as normal.
                    elif acceptor_context_reads:

                        # Determine the acceptor context based on
                        # the reads overlapping the variant.
                        self.debug_msg("ac", variantid)
                        t0 = time.time()
                        acceptor_context = (
                                self.determine_context(acceptor_context_reads,
                                                       vcfvar.get_variant_pos(),
                                                       vcfvar.get_variant_chrom())
                                )
                        self.debug_msg("ac", variantid, t0)

                    # === COMBINED CONTEXT ====================================
                    # Determine the combined variant context based on the
                    # widest window from both the donor and acceptor positions.
                    self.debug_msg("cc", variantid)
                    t0 = time.time()
                    variant_context = (
                            self.determine_largest_context(vcfvar.get_variant_pos(),
                                                           acceptor_context,
                                                           donor_context)
                            )
                    self.debug_msg("cc", variantid, t0)

                    if self.contexts.context_is_in_variant_context(variant_context):
                        self.vaselogger.debug(f"Variant context {variant_context[0]}_{variant_context[1]} overlaps "
                                              f"with an already existing variant context")
                        continue

                    # Obtain all donor reads overlapping the
                    # combined variant context, as well as their mates.
                    self.debug_msg("cdr", variantid)
                    t0 = time.time()
                    variant_context_donor_reads = (
                            self.get_variant_reads(variantid,
                                                   variant_context[0],
                                                   variant_context[2],
                                                   variant_context[3],
                                                   bamfile,
                                                   True,
                                                   varcon_unmapped_d)
                            )
                    self.debug_msg("cdr", variantid, t0)

                    # Obtain all acceptor reads overlapping the
                    # combined variant context, as well as their mates.
                    self.debug_msg("car", variantid)
                    t0 = time.time()
                    variant_context_acceptor_reads = (
                            self.get_variant_reads(variantid,
                                                   variant_context[0],
                                                   variant_context[2],
                                                   variant_context[3],
                                                   acceptorbamfile,
                                                   True,
                                                   varcon_unmapped_a)
                            )
                    self.debug_msg("car", variantid, t0)

                    # If still no acceptor reads were found in the combined
                    # context, set it equal to the dummy set.
                    if not variant_context_acceptor_reads:
                        variant_context_acceptor_reads = None

                    # Add the combined, donor, and acceptor contexts along
                    # with their reads to the current list of contexts.
                    self.contexts.add_variant_context(
                            variantid,
                            sampleid,
                            *variant_context,
                            variant_context_acceptor_reads,
                            variant_context_donor_reads
                            )

                    self.contexts.add_donor_context(
                            variantid,
                            sampleid,
                            *donor_context,
                            donor_context_reads
                            )

                    self.contexts.add_acceptor_context(
                            variantid,
                            sampleid,
                            *acceptor_context,
                            acceptor_context_reads
                            )

                    # Add the read identifiers of reads with
                    # unmapped mates.
                    self.contexts.set_unmapped_acceptor_mate_ids(
                            variantid,
                            varcon_unmapped_a
                            )
                    self.contexts.set_unmapped_donor_mate_ids(
                            variantid,
                            varcon_unmapped_d
                            )
                    self.contexts.set_acceptor_context_unmapped_mate_ids(
                            variantid,
                            acceptor_unmapped
                            )
                    self.contexts.set_donor_context_unmapped_mate_ids(
                            variantid,
                            donor_unmapped
                            )

                except IOError:
                    self.vaselogger.warning("Could not obtain BAM reads from "
                                            f"{bamsamplemap[sampleid]}")
                self.debug_msg("done", variantid, var_t0)
                # End of variant. Proceed to next variant in VCF.

            # End of VCF.
            # Close files and add names of files used.
            bamfile.close()
            donor_vcfs_used.append(vcfsamplemap[sampleid])
            donor_bams_used.append(bamsamplemap[sampleid])
            # Proceed to next sample.

        # End of samples.
        # Write the variant context data and used donor VCFs/BAMs
        # to output files.
        self.vaselogger.info(
                f"Writing combined variant contexts to {varcon_outpath}"
                )
        # Writes the variant contexts to file.
        self.contexts.write_variant_context_file(varcon_outpath)

        self.vaselogger.info(
                "Writing combined variant context statistics to "
                f"{outpath}/varconstats.txt"
                )
        # Writes the variant context statistics file.
        self.contexts.write_variant_context_stats(
                f"{outpath}/varconstats.txt"
                )

        self.vaselogger.info("Writing variant context chrom, start, end, id to a BED file")
        self.write_bed_file(self.contexts.get_variant_contexts(), f"{outpath}/variantcontexts.bed")

        self.vaselogger.info(
                "Write the used donor VCF files per sample to "
                f"{outpath}/donorvcfs.txt"
                )
        self.write_used_donor_files(
                f"{outpath}/donorvcfs.txt",
                vcfsamplemap,
                donor_vcfs_used
                )

        self.vaselogger.info(
                "Write the used donor BAM files per sample to "
                f"{outpath}/donorbams.txt"
                )
        self.write_used_donor_files(
                f"{outpath}/donorbams.txt",
                bamsamplemap,
                donor_bams_used
                )

        # Checks whether the program is running on debug.  If so,
        # write some extra output files.
        if self.vaselogger.getEffectiveLevel() == 10:
            self.write_optional_output_files(outpath, self.contexts)
        acceptorbamfile.close()
        return
# =======

    def build_donor_from_varcon(self, varc_file, bamsamplemap, reference_loc):
        self.contexts = VariantContextFile(varc_file)
        sample_list = []
        for context in self.contexts.variant_contexts.values():
            sample_list.append(context.sample_id)
        sample_list = list(set(sample_list))

        for sampleid in sample_list:
            try:
                bamfile = pysam.AlignmentFile(bamsamplemap[sampleid],
                                              reference_filename=reference_loc)
                self.vaselogger.debug("Opened BAM file "
                                      f"{bamsamplemap[sampleid]}")
            except IOError:
                self.vaselogger.warning("Could not open data files for sample "
                                        f"{sampleid}. Skipping sample.")
                continue
            for context in self.contexts.variant_contexts.values():
                if context.sample_id != sampleid:
                    continue
                context.variant_context_dreads = (
                    self.get_variant_reads(context.context_id,
                                           context.variant_context_chrom,
                                           context.variant_context_start,
                                           context.variant_context_end,
                                           bamfile)
                                           )

    def build_validation_set(self, run_mode,
                             acceptor_bam,
                             fq1_in, fq2_in, fq_out):

        if "F" in run_mode:
            # Set up a set of all acceptor fastq reads to skip.
            skip_list = set(self.contexts.get_all_variant_context_acceptor_read_ids())
            # Set up a list of all donor reads to write.
            add_list = self.contexts.get_all_variant_context_donor_reads()

            # Make the new FastQ files that can be used to run in the
            # NGS_DNA pipeline along real sample data.
            self.vaselogger.info("Start writing the R1 FastQ files.")
            r1fq_starttime = time.time()
            # Build the R1 fastq file.
            self.build_fastq(fq1_in, skip_list, add_list, "F", fq_out)
            self.vaselogger.info("Wrote all R1 FastQ files.")
            self.vaselogger.debug(f"Writing R1 FastQ file(s) took {time.time() - r1fq_starttime} seconds.")

            self.vaselogger.info("Start writing the R2 FastQ files.")
            r2fq_starttime = time.time()
            # Build the R2 fastq file.
            self.build_fastq(fq2_in, skip_list, add_list, "R", fq_out)
            self.vaselogger.info("Wrote all R2 FastQ files.")
            self.vaselogger.debug(f"Writing R2 FastQ file(s) took {time.time() - r2fq_starttime} seconds.")
            return

        elif "D" in run_mode:
            add_list = self.contexts.get_all_variant_context_donor_reads()

            self.vaselogger.info("Only writing donor FastQ files.")

            self.vaselogger.info("Start writing the R1 donor FastQ files.")
            self.build_donor_fq(add_list, fq_out, "F")
            self.vaselogger.info("Finished writing the R1 donor FastQ files.")

            self.vaselogger.info("Start writing the R2 donor FastQ files.")
            self.build_donor_fq(add_list, fq_out, "R")

            self.vaselogger.info("Finished writing the R2 donor FastQ files.")

            self.vaselogger.info("Finished writing donor FastQ files.")
            return

        elif "X" in run_mode:
            return

    # Checks whether a value is in a filter list (array or set)
    def passes_filter(self, val_to_check, filter_to_use, is_exclude_filter=False):
        if filter_to_use is not None:
            if is_exclude_filter:
                if val_to_check in filter_to_use:
                    return False
                return True
            else:
                if val_to_check in filter_to_use:
                    return True
                return False
        return True

    # Returns a list of VcfVariant objects from a VCF file. A filter can be used to select specific variants.
    def get_sample_vcf_variants(self, vcf_fileloc, filterlist=None):
        sample_variant_list = []
        try:
            vcf_file = pysam.VariantFile(vcf_fileloc, "r")
            for vcfvar in vcf_file.fetch():
                if self.passes_filter((vcfvar.chrom, vcfvar.pos), filterlist):
                    varianttype = self.determine_variant_type(vcfvar.ref, vcfvar.alts)
                    sample_variant_list.append(VcfVariant(vcfvar.chrom, vcfvar.pos, vcfvar.ref, vcfvar.alts,
                                                          vcfvar.filter, varianttype))
        except IOError:
            self.vaselogger.warning(f"Could not open VCF file {vcf_fileloc}")
        finally:
            return sample_variant_list

    # Returns whether a variant is a SNP or indel.
    def determine_variant_type(self, vcfvariantref, vcfvariantalts):
        # Determine the maximum reference allele length.
        maxreflength = max([len(x) for x in vcfvariantref.split(",")])
        # Determine the maximum alternative allele length.
        maxaltlength = max([len(x) for x in vcfvariantalts])

        # Check based on the reference and alternative lengths whether
        # the variant is a SNP or indel.
        if maxreflength == 1 and maxaltlength == 1:
            return "snp"
        elif maxreflength > 1 or maxaltlength > 1:
            return "indel"
        return "?"

    # Determines the start and end positions to use for searching reads
    # overlapping with the variant.
    def determine_read_search_window(self, varianttype, vcfvariant):
        if varianttype == "snp":
            return [vcfvariant.get_variant_pos() - 1, vcfvariant.get_variant_pos() + 1]
        elif varianttype == "indel":
            return self.determine_indel_read_range(vcfvariant.get_variant_pos(),
                                                   vcfvariant.get_variant_ref_allele(),
                                                   vcfvariant.get_variant_alt_alleles())
        return [-1, -1]

    # Returns the search start and stop to use for searching BAM reads
    # overlapping with the range of the indel.
    def determine_indel_read_range(self, variantpos, variantref, variantalts):
        searchstart = variantpos
        searchstop = variantpos + max(
                max([len(x) for x in variantref.split(",")]),
                max([len(x) for x in variantalts])
                )
        return [searchstart, searchstop]

    # Returns the BAM reads containing the specific vcf variant as well
    # as their read mate.
    def get_variant_reads(self, contextid, variantchrom,
                          variantstart, variantend,
                          bamfile,
                          write_unm=False, umatelist=None):
        # Obtain all the variant reads overlapping with the variant and
        # their mate reads.
        variantreads = []
        rpnext = {}

        for vread in bamfile.fetch(variantchrom, variantstart, variantend):
            variantreads.append(DonorBamRead(
                    vread.query_name,
                    self.get_read_pair_num(vread),
                    vread.reference_name,
                    vread.reference_start,
                    vread.infer_read_length(),
                    vread.get_forward_sequence(),
                    "".join([chr((x + 33))
                             for x in vread.get_forward_qualities()]),
                    vread.mapping_quality
                    ))
            rpnext[vread.query_name] = [vread.next_reference_name, vread.next_reference_start, vread.query_name]

        # Obtain the read mate (this must be done after the initial fetch not during!)
        for read1 in rpnext.values():
            materead = self.fetch_mate_read(*read1, bamfile)
            if materead is not None:
                variantreads.append(materead)
            else:
                if write_unm:
                    self.vaselogger.debug("Could not find mate for "
                                          f"{read1[2]} ; "
                                          "mate is likely unmapped.")
                    umatelist.append(read1[2])

        # Make sure the list only contains each BAM read once (if a read
        # and mate both overlap with a variant, they have been added
        # twice to the list).
        uniq_variantreads = []
        checklist = []
        for fetched in variantreads:
            id_pair = (fetched.get_bam_read_id(), fetched.get_bam_read_pair_number())
            if id_pair not in checklist:
                uniq_variantreads.append(fetched)
                checklist.append(id_pair)
        variantreads = uniq_variantreads

        # Filter to keep only read pairs.
        variantreads = self.filter_variant_reads(variantreads)
        self.vaselogger.debug(f"Found a total of {len(variantreads)} "
                              "BAM reads.")
        return variantreads

    # Returns the mate read using the fetch() method
    def fetch_mate_read(self, rnext, pnext, readid, bamfile):
        for bamread in bamfile.fetch(rnext, pnext, pnext + 1):
            if bamread.query_name == readid and bamread.reference_start == pnext:
                return DonorBamRead(bamread.query_name, self.get_read_pair_num(bamread), bamread.reference_name,
                                    bamread.reference_start, bamread.infer_read_length(),
                                    bamread.get_forward_sequence(),
                                    "".join([chr((x + 33)) for x in bamread.get_forward_qualities()]),
                                    bamread.mapping_quality)
        return None

    # Filters the donor reads to keep only reads that occur twice.
    def filter_variant_reads(self, bamreads):
        return [bread for bread in bamreads
                if (self.read_occurence(bread.get_bam_read_id(), bamreads) == 2)]

    # Returns the number of occurences of a certain read in the list of
    # BAM reads (should be two ideally).
    def read_occurence(self, readid, readlist):
        return sum([bamread.get_bam_read_id() == readid for bamread in readlist])

    # Determines the start and stops of the variant context (please see
    # the documentation for more information).
    def determine_context(self, contextreads, contextorigin, contextchr):
        # Check whether there are reads to determine the context for.
        if not contextreads:
            return []

        # Get read start and stop position, while filtering out
        # read mates that map to different chr.
        starts = [conread.get_bam_read_ref_pos()
                  for conread in contextreads
                  if conread.get_bam_read_chrom() == contextchr]

        stops = [conread.get_bam_read_ref_end()
                 for conread in contextreads
                 if conread.get_bam_read_chrom() == contextchr]

        # Number of mates filtered for mapping to different chr.
        num_diff_chr = len(contextreads) - len(starts)

        # Filter outlier read start/stops.
        filtered_starts = self.filter_outliers(starts)
        filtered_stops = self.filter_outliers(stops)

        # Number of filtered outlier positions.
        num_filtered = (len(starts)
                        + len(stops)
                        - len(filtered_starts)
                        - len(filtered_stops))

        self.vaselogger.debug(f"{num_diff_chr} read mate(s) filtered due to "
                              "alignment to different reference sequence.")
        self.vaselogger.debug(f"{num_filtered} outlier read position(s) filtered.")

        # Set variant context as chr, min start, max end.
        contextstart = min(filtered_starts)
        contextend = max(filtered_stops)
        self.vaselogger.debug("Context is "
                              + str(contextchr) + ", "
                              + str(contextstart) + ", "
                              + str(contextend))
        return [contextchr, contextorigin, contextstart, contextend]

    # Filters outliers from a list of start/stop positions using
    # Tukey's Fences method.
    def filter_outliers(self, pos_list, k=3):
        # First and third quartile values of the positions.
        q1 = np.percentile(pos_list, 25)
        q3 = np.percentile(pos_list, 75)
        # Interquartile range.
        iq = q3 - q1
        # Only include positions that fall within the range of
        # (q1 to q3) +/- k*iq.
        filtered = [x for x in pos_list
                    if ((q1 - (k * iq)) <= x <= (q3 + (k * iq)))]
        return filtered

    # Determines the size of the variant context based on both the
    # acceptor and donor reads.
    def determine_largest_context(self, contextorigin, acceptor_context,
                                  donor_context):
        largest_context = [donor_context[0]]
        largest_context.append(contextorigin)
        # Determine and save the smallest context start.
        largest_context.append(min(acceptor_context[2], donor_context[2]))
        # Determine and save the largest context end.
        largest_context.append(max(acceptor_context[3], donor_context[3]))
        return largest_context

    # ===METHODS TO PRODUCE THE VALIDATION FASTQ FILES=========================
    # Will build the R1/R2 VaSe fastq files.
    def build_fastq(self, acceptorfq_filepaths, acceptorreads_toskip,
                    donor_context_reads, forward_or_reverse, vasefq_outpath):
        writedonor = False

        # Iterate over the R1/R2 fastq in files to use as templates for the
        for x in range(0, len(acceptorfq_filepaths)):
            if x == len(acceptorfq_filepaths)-1:
                writedonor = True
                self.vaselogger.debug("Donor reads will be added the current "
                                      "VaSe fastQ out file.")

            # Write the new VaSe FastQ file.
            vasefq_outname = self.set_fastq_out_path(vasefq_outpath,
                                                     forward_or_reverse,
                                                     x + 1)
            self.vaselogger.debug(f"Set FastQ output path to: {vasefq_outname}")
            self.write_vase_fastq(acceptorfq_filepaths[x], vasefq_outname,
                                  acceptorreads_toskip, donor_context_reads,
                                  forward_or_reverse, writedonor)

    # Builds a new FastQ file to be used for validation.
    def write_vase_fastq(self, acceptor_infq, fastq_outpath,
                         acceptorreads_toskip, donorbamreaddata,
                         fr, writedonordata=False):
        try:
            fqgz_outfile = io.BufferedWriter(open(fastq_outpath, "wb"))
            self.vaselogger.debug(f"Opened template FastQ: {acceptor_infq}")

            # Open the template fastq and write filtered data to a new
            # fastq.gz file.
            fqgz_infile = io.BufferedReader(gzip.open(acceptor_infq, "rb"))
            for fileline in fqgz_infile:

                # Check if we are located at a read identifier.
                if fileline.startswith(b"@"):
                    if fileline.decode("utf-8").split()[0][1:] not in acceptorreads_toskip:
                        fqgz_outfile.write(fileline)
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
            fqgz_infile.close()

            # Add the patient BAM reads containing a VCF variant to the
            # new FastQ file.
            if writedonordata:
                donorbamreaddata.sort(key=lambda dbr: dbr[0], reverse=False)
                for bamread in donorbamreaddata:
                    # Check if the BAM read is R1 or R2.
                    if bamread[1] == fr:
                        fqlines = ("@" + str(bamread[0]) + "\n"
                                   + str(bamread[2]) + "\n"
                                   + "+\n"
                                   + str(bamread[3]) + "\n")
                        fqgz_outfile.write(fqlines.encode("utf-8"))
            fqgz_outfile.flush()
            fqgz_outfile.close()

        except IOError as ioe:
            if ioe.filename == acceptor_infq:
                self.vaselogger.critical("The supplied template FastQ file "
                                         "could not be found.")
            if ioe.filename == fastq_outpath:
                self.vaselogger.critical("A FastQ file could not be written "
                                         "to the provided output location.")
            exit()

    def build_donor_fq(self, donorbamreaddata, fastq_outpath, fr):
        # Write the new VaSe FastQ file.
        vasefq_outname = self.set_fastq_out_path(fastq_outpath, fr, 1)
        fqgz_outfile = io.BufferedWriter(open(vasefq_outname, "wb"))

        donorbamreaddata.sort(key=lambda dbr: dbr[0], reverse=False)
        for bamread in donorbamreaddata:
            # Check if the BAM read is R1 or R2.
            if bamread[1] == fr:
                fqlines = ("@" + str(bamread[0]) + "\n"
                           + str(bamread[2]) + "\n"
                           + "+\n"
                           + str(bamread[3]) + "\n")
                fqgz_outfile.write(fqlines.encode("utf-8"))
        fqgz_outfile.flush()
        fqgz_outfile.close()

    # Checks if a read is read 1 (R1) or read 2 (R2).
    def is_required_read(self, bamread, fr):
        if fr == "F":
            return bamread.is_read1()
        return bamread.is_read2()

    # Returns the name for the fastq out file.
    def set_fastq_out_path(self, outpath, fr, lnum):
        if fr == "1":
            return f"{outpath}_{datetime.now().date()}_L{lnum}_R1.fastq"
        return f"{outpath}_{datetime.now().date()}_L{lnum}_R2.fastq"

    # ===METHODS TO OBTAIN SOME DATA OF THE VASEBUILDER OBJECT=================
    # Returns the identifier of the current VaSeBuilder object.
    def get_creation_id(self):
        return self.creation_id

    # Returns the date the current VaSeBuilder object has been made.
    def get_creation_date(self):
        return self.creation_date

    # Returns the time the current VaSeBuilder object has been made.
    def get_creation_time(self):
        return self.creation_time

    # Returns all variant contexts.
    def get_variant_contexts(self):
        return self.contexts

    # Returns a specified acceptor context.
    def get_acceptor_context(self, contextid):
        return self.contexts.get_acceptor_context(contextid)

    # Returns a specified donor context.
    def get_donor_context(self, contextid):
        return self.contexts.get_donor_context(contextid)

    # Returns the context start and stop for a specified VCF variant.
    def get_variant_context(self, contextid):
        return self.contexts.get_variant_context(contextid)

    # Returns an identifier for a VCF variant.
    # If the identifier is "." then one will be constructed as
    # "chrom_pos".
    def get_vcf_variant_id(self, vcfvariant):
        return f"{vcfvariant.chrom}_{vcfvariant.pos}"

    # Returns whether the read is the first or second read in a pair.
    def get_read_pair_num(self, pysam_bamread):
        if (pysam_bamread.is_read1):
            return "1"
        return "2"

    # Returns the chromosome prefix to use for a BAM/CRAM file ('chr' or '')
    def get_searchchrom_prefix(self, bamfile):
        if bamfile.header["SQ"][0]["SN"].startswith('chr'):
            return "chr"
        elif bamfile.header["SQ"][0]["SN"].startswith("CHR"):
            return "CHR"
        return ""

    # ===METHODS TO WRITE OUTPUT FILES=========================================
    # Writes the used donor vcf files to a file
    def write_used_donor_files(self, outfileloc, filesamplemap,
                               used_donor_files):
        try:
            with open(outfileloc, "w") as outfile:
                outfile.write("#SampleId\tDonorFile\n")
                for sampleid, samplefile in filesamplemap.items():
                    if (samplefile in used_donor_files):
                        outfile.write(f"{sampleid}\t{samplefile}\n")
        except IOError as ioe:
            self.vaselogger.critical("Could not write used donor files to "
                                     f"{outfileloc}")

    # Writes the optional output files (when logger is set to DEBUG log
    # level).
    def write_optional_output_files(self, outpath, contextfile):
        # Write the optional acceptor context files; acceptor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaselogger.debug("Writing acceptor contexts to "
                              f"{outpath}/acceptorcontexts.txt")
        contextfile.write_acceptor_context_file(f"{outpath}/acceptorcontexts.txt")

        self.vaselogger.debug("Writing acceptor context statistics to "
                              f"{outpath}/acceptorcontextstats.txt")
        contextfile.write_acceptor_context_stats(
                              f"{outpath}/acceptorcontextstats.txt"
                              )

        self.vaselogger.debug("Writing acceptor context read identifiers with "
                              "unmapped mates to"
                              f"mates to {outpath}/acceptor_unmapped.txt")
        contextfile.write_acceptor_unmapped_mates(
                              f"{outpath}/acceptor_unmapped.txt"
                              )

        self.vaselogger.debug("Writing left and right most read positions of "
                              "each acceptor context to "
                              f"{outpath}/acceptor_positions.txt")
        contextfile.write_acceptor_left_right_positions(
                              f"{outpath}/acceptor_positions.txt"
                              )

        # Write the optional donor context files; donor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaselogger.debug("Writing donor contexts to "
                              f"{outpath}/donorcontexts.txt")
        contextfile.write_donor_context_file(f"{outpath}/donorcontexts.txt")

        self.vaselogger.debug("Writing donor context statistics to "
                              f"{outpath}/donorcontextstats.txt")
        contextfile.write_donor_context_stats(f"{outpath}/donorcontextstats.txt")

        self.vaselogger.debug("Writing donor context read identifiers "
                              "with unmapped mates to "
                              f"{outpath}/donor_unmapped.txt")
        contextfile.write_donor_unmapped_mates(f"{outpath}/donor_unmapped.txt")

        self.vaselogger.debug("Writing left and right most read positions "
                              "of each donor context to "
                              f"{outpath}/donor_positions.txt")
        contextfile.write_donor_left_right_positions(
                              f"{outpath}/donor_positions.txt"
                              )

        # Write the optional variant context files; acceptor & donor
        # unmapped mates and left/right positions.
        self.vaselogger.debug("Writing variant context acceptor read "
                              "identifiers with unmapped mates to "
                              f"{outpath}/varcon_unmapped_acceptor.txt")
        contextfile.write_reads_with_unmapped_mate(
                              "acceptor",
                              f"{outpath}/varcon_unmapped_acceptor.txt"
                              )

        self.vaselogger.debug("Writing variant context donor read identifiers "
                              "with unmapped mates to "
                              f"{outpath}/varcon_unmapped_donor.txt")
        contextfile.write_reads_with_unmapped_mate(
                              "donor",
                              f"{outpath}/varcon_unmapped_donor.txt"
                              )

        self.vaselogger.debug("Writing variant context left and right most "
                              "read positions of acceptor reads to "
                              f"{outpath}/varcon_positions_acceptor.txt")
        contextfile.write_left_right_positions(
                              "acceptor",
                              f"{outpath}/varcon_positions_acceptor.txt"
                              )

        self.vaselogger.debug("Writing variant context left and right most "
                              "read positions of donor reads to "
                              f"{outpath}/varcon_positions_donor.txt")
        contextfile.write_left_right_positions(
                              "donor",
                              f"{outpath}/varcon_positions_donor.txt"
                              )

    # Writes a BED file for the variant context data
    def write_bed_file(self, variantcontextdata, bedoutloc):
        try:
            with open(bedoutloc, "w") as bedoutfile:
                for varcon in variantcontextdata:
                    bedoutfile.write(f"{varcon.get_variant_context_chrom()}\t{varcon.get_variant_context_start()}\t"
                                     f"{varcon.get_variant_context_end()}\t{varcon.get_variant_context_id()}\n")
        except IOError:
            self.vaselogger.warning(f"Could not write variant context data to BED file: {bedoutloc}")
