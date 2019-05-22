#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import statistics
import pysam

# Import VaSe specific classes.
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile


class VaSeBuilder:
    # Constructor that saves the identifier, date and time of the current
    def __init__(self, vaseId):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.creation_id = str(vaseId)
        self.creation_date = datetime.now().date()
        self.creation_time = datetime.now().time()

        # Create the bookkeeping variables used for saving the variant
        # contexts.
        # VariantContextFile that saves the acceptor, donor and variant
        # contexts and their associated data.
        self.contexts = VariantContextFile()
        self.vaselogger.info("VaSeBuilder: "
                             + str(self.creation_id) + " ; "
                             + str(self.creation_date) + " ; "
                             + str(self.creation_time))

    # ===METHODS TO PRODUCE VARIANT CONTEXTS===================================
    # Creates the new FastQ validation dataset by replacing NIST reads
    # containing a VCF variant with BAM reads from patients.  Returns
    # true at the end to indicate the process is done.
    def build_validation_set(self, vcfBamLinkMap,
                             vcfsamplemap, bamsamplemap,
                             acceptorbamloc,
                             fastq_fpath, fastq_rpath,
                             outpath,
                             fastq_outpath,
                             varcon_outpath):
        self.vaselogger.info("Start building the validation set")
        donor_vcfs_used, donor_bams_used = [], []

        try:
            acceptorbamfile = pysam.AlignmentFile(acceptorbamloc, 'rb')

            # Iterate over the samples to use for building the validation set.
            for sampleid in vcfsamplemap:
                self.vaselogger.debug("Processing data for sample "
                                      + str(sampleid))

                # Check in the VCF BAM link map that there is a BAM file
                # for the sample as well.
                try:
                    vcffile = pysam.VariantFile(vcfsamplemap[sampleid], 'r')
                    self.vaselogger.debug("Opened VCF file "
                                          + str(vcfsamplemap[sampleid]))
                    bamfile = pysam.AlignmentFile(bamsamplemap[sampleid], 'rb')
                    self.vaselogger.debug("Opened BAM file "
                                          + str(bamsamplemap[sampleid]))

                    # Loop over the variants in the VCF file. Prior to
                    # identifying the BAM reads, it is first checked
                    # whether the variant is in a previously established
                    for vcfvar in vcffile.fetch():
                        variantid = self.get_vcf_variant_id(vcfvar)
                        acceptor_unmapped = []
                        donor_unmapped = []
                        varcon_unmapped_a = []
                        varcon_unmapped_d = []
                        self.vaselogger.debug(
                                "Searching BAM reads for variant "
                                + str(variantid)
                                )

                        # Determine the search window before gathering
                        # reads overlapping with the VCF variant.
                        varianttype = self.determine_variant_type(vcfvar.ref,
                                                                  vcfvar.alts)
                        searchwindow = self.determine_read_search_window(
                                varianttype,
                                vcfvar
                                )

                        # Get the BAM reads for the variant and
                        # determine the variant context.
                        if (not self.contexts.variant_is_in_context(
                                varianttype,
                                vcfvar.chrom,
                                searchwindow[0],
                                searchwindow[1]
                                )):
                            try:
                                # Gather acceptor reads and their mates
                                # overlapping with the variant and determine
                                # the acceptor context.
                                self.vaselogger.debug(
                                        "Determine acceptor BAM reads for "
                                        "variant " + str(variantid)
                                        )
                                # Obtain all acceptor BAM reads containing the
                                # VCF variant and their read mate.
                                acceptor_context_reads = self.get_variant_reads(
                                        variantid, vcfvar.chrom,
                                        searchwindow[0], searchwindow[1],
                                        acceptorbamfile, True,
                                        acceptor_unmapped
                                        )
                                self.vaselogger.debug("Determine acceptor "
                                                      "context for variant "
                                                      + str(variantid))
                                # Determine the acceptor variant context based
                                # on the reads overlapping the variant.
                                acceptor_context = self.determine_context(
                                        acceptor_context_reads,
                                        vcfvar.pos
                                        )

                                # Gather donor reads and their mates
                                # overlapping with the variant and determine
                                # the donor context.
                                self.vaselogger.debug("Search donor BAM reads "
                                                      "for variant "
                                                      + str(variantid))
                                # Obtain all donot BAM reads containing the VCF
                                # variant and their read mate.
                                donor_context_reads = self.get_variant_reads(
                                        variantid, vcfvar.chrom,
                                        searchwindow[0], searchwindow[1],
                                        bamfile, True,
                                        donor_unmapped
                                        )
                                self.vaselogger.debug("Determine donor "
                                                      "context for variant "
                                                      + str(variantid))
                                # Determine the donor variant context based on
                                # the reads overlapping the variant.
                                donor_context = self.determine_context(
                                        donor_context_reads,
                                        vcfvar.pos
                                        )

                                # Determine the ultimate variant context and
                                # obtain the overlapping acceptor and donor
                                # reads.
                                variant_context = self.determine_largest_context(
                                        vcfvar.pos,
                                        acceptor_context,
                                        donor_context
                                        )
                                # Obtain all acceptor reads overlapping with
                                # the combined variant context and their mates.
                                variant_context_acceptor_reads = self.get_variant_reads(
                                        variantid,
                                        variant_context[0],
                                        variant_context[2],
                                        variant_context[3],
                                        acceptorbamfile, True,
                                        varcon_unmapped_a
                                        )
                                # Obtain all donor reads overlapping with the
                                # combined variant context and their mates.
                                variant_context_donor_reads = self.get_variant_reads(
                                        variantid,
                                        variant_context[0],
                                        variant_context[2],
                                        variant_context[3],
                                        bamfile, True,
                                        varcon_unmapped_d
                                        )

                                # Check whether reads were found in both
                                # acceptor and donor.  Only then save the
                                # results.
                                if ((len(donor_context_reads) > 0)
                                   and (len(acceptor_context_reads) > 0)):
                                    self.contexts.add_variant_context(
                                            variantid,
                                            sampleid,
                                            variant_context[0],
                                            variant_context[1],
                                            variant_context[2],
                                            variant_context[3],
                                            variant_context_acceptor_reads,
                                            variant_context_donor_reads
                                            )
                                    self.contexts.add_acceptor_context(
                                            variantid,
                                            sampleid,
                                            acceptor_context[0],
                                            acceptor_context[1],
                                            acceptor_context[2],
                                            acceptor_context[3],
                                            acceptor_context_reads
                                            )
                                    self.contexts.add_donor_context(
                                            variantid,
                                            sampleid,
                                            donor_context[0],
                                            donor_context[1],
                                            donor_context[2],
                                            donor_context[3],
                                            donor_context_reads
                                            )

                                    # Add the read identifiers of reads with
                                    # an unmapped mate.
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
                                else:
                                    self.vaselogger.debug(
                                            "No donor and/or acceptor BAM "
                                            "reads found for variant "
                                            + str(variantid)
                                            )

                            except IOError as ioe:
                                self.vaselogger.warning(
                                        "Could not obtain BAM reads from "
                                        + str(bamsamplemap[sampleid])
                                        )
                        else:
                            self.vaselogger.debug(
                                    "VCF variant "
                                    + str(variantid)
                                    + " is located in an already existing "
                                    "variant context"
                                    )
                    bamfile.close()
                    vcffile.close()

                    donor_vcfs_used.append(vcfsamplemap[sampleid])
                    donor_bams_used.append(bamsamplemap[sampleid])
                except IOError as ioe:
                    self.vaselogger.warning("Could not establish data for "
                                            + str(sampleid))

            # Write the variant context data and used donor VCFs/BAMs
            # to output files.
            self.vaselogger.info(
                    "Writing combined variant contexts to "
                    + str(varcon_outpath)
                    )
            # Writes the variant contexts to file.
            self.contexts.write_variant_context_file(varcon_outpath)

            self.vaselogger.info(
                    "Writing combined variant context statistics to "
                    + str(outpath) + "varconstats.txt"
                    )
            # Writes the variant context statistics file.
            self.contexts.write_variant_context_stats(
                outpath + "/varconstats.txt"
                    )

            self.vaselogger.info(
                    "Write the used donor VCF files per sample to "
                    + str(outpath) + "/donorvcfs.txt"
                    )
            self.write_used_donor_files(
                outpath + "/donorvcfs.txt",
                    vcfsamplemap,
                    donor_vcfs_used
                    )

            self.vaselogger.info(
                    "Write the used donor BAM files per sample to "
                    + str(outpath) + "/donorbams.txt"
                    )
            self.write_used_donor_files(
                outpath + "/donorbams.txt",
                    bamsamplemap,
                    donor_bams_used
                    )

            # Checks whether the program is running on debug.  If so,
            # write some extra output files.
            if self.vaselogger.getEffectiveLevel() == 10:
                self.write_optional_output_files(outpath, self.contexts)

            # Obtain a list of acceptor reads to skip when iterating
            # over the acceptor FastQ.
            # Set up a list of all acceptor reads to skip.
            acceptor_reads_to_skip = list(set(
                    self.contexts.get_all_variant_context_acceptor_read_ids()
                    ))
            # Sets up a list.
            donorreads = self.contexts.get_all_variant_context_donor_reads()

            # Make the new FastQ files that can be used to run in the
            # NGS_DNA pipeline along real sample data.
            self.vaselogger.info("Start writing the R1 FastQ files")
            # Build the R1 fastq file.
            self.build_fastq(fastq_fpath, acceptor_reads_to_skip,
                             donorreads, 'F', fastq_outpath)
            self.vaselogger.info("Wrote all R1 FastQ files")

            self.vaselogger.info("Start writing the R2 FastQ files")
            # Build the R2 fastq file.
            self.build_fastq(fastq_rpath, acceptor_reads_to_skip,
                             donorreads, 'R', fastq_outpath)
            self.vaselogger.info("Wrote all R2 FastQ files")

            self.vaselogger.info("Finished building the validation set")
            acceptorbamfile.close()

        except IOError as ioe:
            self.vaselogger.critical("Could not open acceptor BAM file")
            exit()

    # Returns whether a variant is a SNP or indel.
    def determine_variant_type(self, vcfvariantref, vcfvariantalts):
        # Determine the maximum reference allele length.
        maxreflength = max([len(x) for x in vcfvariantref.split(',')])
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
            return [vcfvariant.pos - 1, vcfvariant.pos + 1]
        elif varianttype == "indel":
            return self.determine_indel_read_range(vcfvariant.pos,
                                                   vcfvariant.ref,
                                                   vcfvariant.alts)
        return [-1, -1]

    # Returns the search start and stop to use for searching BAM reads
    # overlapping with the range of the indel.
    def determine_indel_read_range(self, variantpos, variantref, variantalts):
        searchstart = variantpos
        searchstop = variantpos + max(
                max([len(x) for x in variantref.split(',')]),
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

        for vread in bamfile.fetch(variantchrom, variantstart, variantend):
            variantreads.append(DonorBamRead(
                    vread.query_name,
                    self.get_read_pair_num(vread),
                    vread.reference_name,
                    vread.reference_start,
                    vread.infer_read_length(),
                    vread.get_forward_sequence(),
                    ''.join([chr((x + 33))
                             for x in vread.get_forward_qualities()]),
                    vread.mapping_quality
                    ))

            # Try to obtain the reads mate as well.
            try:
                vmread = bamfile.mate(vread)
                variantreads.append(DonorBamRead(
                        vmread.query_name,
                        self.get_read_pair_num(vmread),
                        vmread.reference_name,
                        vmread.reference_start,
                        vmread.infer_read_length(),
                        vmread.get_forward_sequence(),
                        ''.join([chr((x+33))
                                 for x in vmread.get_forward_qualities()]),
                        vmread.mapping_quality
                        ))

            except ValueError as pve:
                self.vaselogger.debug("Could not find mate for "
                                      + vread.query_name
                                      + " ; mate is likely unmapped.")
                if write_unm:
                    umatelist.append(vread.query_name)

        # Make sure the list only contains each BAM read once (if a read
        # and mate both overlap with a variant, they have been added
        # twice to the list).
        variantreads = list(set(variantreads))
        variantreads = self.filter_variant_reads(variantreads)
        self.vaselogger.debug("Found a total of "
                              + str(len(variantreads))
                              + " BAM reads.")
        return variantreads

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
    def determine_context(self, contextreads, contextorigin):
        # Check whether there are reads to determine the context for.
        if len(contextreads) > 0:
            contextchrom = contextreads[0].get_bam_read_chrom()
            contextstart = min([conread.get_bam_read_ref_pos()
                                for conread in contextreads])
            contextend = max([conread.get_bam_read_ref_end()
                              for conread in contextreads])
            self.vaselogger.debug("Context is "
                                  + str(contextchrom) + ", "
                                  + str(contextstart) + ", "
                                  + str(contextend))
            return [contextchrom, contextorigin, contextstart, contextend]
        return []

    # Determines the size of the variant context based on both the
    # acceptor and donor reads.
    def determine_largest_context(self, contextorigin, acceptor_context,
                                  donor_context):
        largest_context = [acceptor_context[0]]
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
            self.vaselogger.debug("Set FastQ output path to: "
                                  + str(vasefq_outname))
            self.write_vase_fastq(acceptorfq_filepaths[x], vasefq_outname,
                                  acceptorreads_toskip, donor_context_reads,
                                  forward_or_reverse, writedonor)

    # Builds a new FastQ file to be used for validation.
    def write_vase_fastq(self, acceptor_infq, fastq_outpath,
                         acceptorreads_toskip, donorbamreaddata,
                         fr, writedonordata=False):
        try:
            fqgz_outfile = io.BufferedWriter(gzip.open(fastq_outpath, 'wb'))
            self.vaselogger.debug("Opened template FastQ: "
                                  + str(acceptor_infq))

            # Open the template fastq and write filtered data to a new
            # fastq.gz file.
            fqgz_infile = io.BufferedReader(gzip.open(acceptor_infq, 'rb'))
            for fileline in fqgz_infile:

                # Check if we are located at a read identifier.
                if fileline.startswith(b"@"):
                    if fileline.decode("utf-8").strip()[1:] not in acceptorreads_toskip:
                        fqgz_outfile.write(fileline)
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
            fqgz_infile.close()

            # Add the patient BAM reads containing a VCF variant to the
            # new FastQ file.
            if writedonordata:
                donorbamreaddata.sort(key=lambda x: x.get_bam_read_id(),
                                      reverse=False)
                for bamread in donorbamreaddata:
                    # Check if the BAM read is R1 or R2.
                    if self.is_required_read(bamread, fr):
                        fqgz_outfile.write(bamread.get_as_fastq_seq().encode("utf-8"))
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

    # Checks if a read is read 1 (R1) or read 2 (R2).
    def is_required_read(self, bamread, fr):
        if fr == "F":
            return bamread.is_read1()
        return bamread.is_read2()

    # Returns the name for the fastq out file.
    def set_fastq_out_path(self, outpath, fr, lnum):
        if fr == "F":
            return (str(outpath)
                    + "_" + str(datetime.now().date())
                    + "_L" + str(lnum)
                    + "_R1.fastq.gz")
        return (str(outpath)
                + "_" + str(datetime.now().date())
                + "_L" + str(lnum)
                + "_R2.fastq.gz")

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
    # If the identifier is '.' then one will be constructed as
    # 'chrom_pos'.
    def get_vcf_variant_id(self, vcfvariant):
        return str(vcfvariant.chrom) + "_" + str(vcfvariant.pos)

    # Returns whether the read is the first or second read in a pair.
    def get_read_pair_num(self, pysam_bamread):
        if (pysam_bamread.is_read1):
            return '1'
        return '2'

    # ===METHODS TO WRITE OUTPUT FILES=========================================
    # Writes the used donor vcf files to a file
    def write_used_donor_files(self, outfileloc, filesamplemap,
                               used_donor_files):
        try:
            with open(outfileloc, 'w') as outfile:
                outfile.write("#SampleId\tDonorFile")
                for sampleid, samplefile in filesamplemap.items():
                    if (samplefile in used_donor_files):
                        outfile.write(sampleid + "\t" + str(samplefile) + "\n")
        except IOError as ioe:
            self.vaselogger.critical("Could not write used donor files to "
                                     + str(outfileloc))

    # Writes the optional output files (when logger is set to DEBUG log
    # level).
    def write_optional_output_files(self, outpath, contextfile):
        # Write the optional acceptor context files; acceptor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaselogger.debug(
                "Writing acceptor contexts to "
                + str(outpath) + "/acceptorcontexs.txt"
                )
        contextfile.write_acceptor_context_file(
            outpath + "/acceptorcontexts.txt"
                )
        self.vaselogger.debug(
                "Writing acceptor context statistics to "
                + str(outpath) + "/acceptorcontextstats.txt"
                )
        contextfile.write_acceptor_context_stats(
            outpath + "/acceptorcontextstats.txt"
                )
        self.vaselogger.debug(
                "Writing acceptor context read identifiers "
                "with unmapped mates to "
                + str(outpath) + "/acceptor_unmapped.txt"
                )
        contextfile.write_acceptor_unmapped_mates(
            outpath + "/acceptor_unmapped.txt"
                )
        self.vaselogger.debug(
                "Writing left and right most read positions "
                "of each acceptor context to "
                + str(outpath) + "/acceptor_positions.txt"
                )
        contextfile.write_acceptor_left_right_positions(
            outpath + "/acceptor_positions.txt"
                )

        # Write the optional donor context files; donor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaselogger.debug(
                "Writing donor contexts to "
                + str(outpath) + "/donorcontexs.txt"
                )
        contextfile.write_donor_context_file(
            outpath + "/donorcontexts.txt"
                )
        self.vaselogger.debug(
                "Writing donor context statistics to "
                + str(outpath) + "/'donorcontextstats.txt"
                )
        contextfile.write_donor_context_stats(
            outpath + "/donorcontextstats.txt"
                )
        self.vaselogger.debug(
                "Writing donor context read identifiers "
                "with unmapped mates to "
                + str(outpath) + "/donor_unmapped.txt"
                )
        contextfile.write_donor_unmapped_mates(
            outpath + "/donor_unmapped.txt"
                )
        self.vaselogger.debug(
                "Writing left and right most read positions "
                "of each donor context to "
                + str(outpath) + "/donor_positions.txt"
                )
        contextfile.write_donor_left_right_positions(
            outpath + "/donor_positions.txt"
                )

        # Write the optional variant context files; acceptor & donor
        # unmapped mates and left/right positions.
        self.vaselogger.debug(
                "Writing variant context acceptor read identifiers "
                "with unmapped mates to "
                + str(outpath) + "/varcon_unmapped_acceptor.txt"
                )
        contextfile.write_reads_with_unmapped_mate(
                'acceptor',
            outpath + "/varcon_unmapped_acceptor.txt"
                )
        self.vaselogger.debug(
                "Writing variant context donor read identifiers "
                "with unmapped mates to "
                + str(outpath) + "/varcon_unmapped_donor.txt"
                )
        contextfile.write_reads_with_unmapped_mate(
                'donor',
            outpath + "/varcon_unmapped_donor.txt"
                )

        self.vaselogger.debug(
                "Writing variant context left and right most read positions "
                "of acceptor reads to "
                + str(outpath) + "/varcon_positions_acceptor.txt"
                )
        contextfile.write_left_right_positions(
                'acceptor',
            outpath + "/varcon_positions_acceptor.txt"
                )
        self.vaselogger.debug(
                "Writing variant context left and right most read positions "
                "of donor reads to "
                + str(outpath) + "/varcon_positions_donor.txt"
                )
        contextfile.write_left_right_positions(
                'donor',
            outpath + "/varcon_positions_donor.txt"
                )
