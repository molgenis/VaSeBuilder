#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import numpy as np
import pysam
import time
import random

# Import VaSe specific classes.
from VcfBamScanner import VcfBamScanner
from DonorBamRead import DonorBamRead
from VariantContextFile import VariantContextFile
from VcfVariant import VcfVariant
from VariantContext import VariantContext
from OverlapContext import OverlapContext


class VaSeBuilder:
    """Creates the variant contexts, builds validation sets and writes the output files.

    Attributes
    ----------
    vb_scanner : VcfBamScanner
        Scans and extracts info from variant and alignment files
    creation_id : str
        Unique (uuid) identifier to identify VaSeBuilder runs
    creation_time : str
        The date and time of creation to identify VaSeBuilder runs
    vaselogger : Logger
        VaSeBuilder logger to log VaSeBuilder activity
    """
    # Constructor that saves the identifier, date, and time of the current run.
    def __init__(self, vaseid):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vb_scanner = VcfBamScanner()
        self.creation_id = str(vaseid)
        self.creation_time = datetime.now()
        self.vaselogger.info(
            f"VaSeBuilder: {self.creation_id} ; {self.creation_time}"
        )

        # VariantContextFile that saves the acceptor, donor, and variant contexts with their associated data.
        self.contexts = VariantContextFile()

        # Dictionary used for debug messages.
        self.debug_dict = {"vw": "Establishing search window",
                           "dr": "Gathering donor variant reads",
                           "dc": "Determining donor context",
                           "ar": "Gathering acceptor variant reads",
                           "ac": "Determining acceptor context",
                           "cc": "Determining combined context",
                           "cdr": "Gathering combined context donor reads",
                           "car": "Gathering combined context acceptor reads",
                           "done": "Variant complete. Processing"}

        self.cigar_tuple_table = {"M": 0,
                                  "I": 1,
                                  "D": 2,
                                  "N": 3,
                                  "S": 4,
                                  "H": 5,
                                  "P": 6,
                                  "-": 7,
                                  "X": 8,
                                  "B": 9}

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
    def passes_filter(self, val_to_check, filter_to_use, is_exclude_filter=False):
        """Checks whether a value is in a filter list.

        A value can be checked against either an inclusion or exclusion filter. An inclusion filter should be the values
        that should be included and used, an exclusion filter for values to be excluded and not used.

        Parameters
        ----------
        val_to_check : str
            Value to check against filter
        filter_to_use : list of str
            Values to use as filter
        is_exclude_filter : bool
            Is the filter an exclusion filter (false for inclusion filter)

        Returns
        -------
        bool
            True if value in an inclusive filter or not in an exclusive filter, False otherwise
        """
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
        """Reads and returns read variants from a variant file.

        Parameters
        ----------
        vcf_fileloc : str
            Path to variant file to read
        filterlist : list of tuple
            Variants to include

        Returns
        -------
        sample_variant_list : list of pysam.VariantRecord
            Read variants fro the variant file
        """
        sample_variant_list = []
        try:
            vcf_file = pysam.VariantFile(vcf_fileloc, "r")
            for vcfvar in vcf_file.fetch():
                if self.passes_filter((vcfvar.chrom, vcfvar.pos), filterlist):
                    sample_variant_list.append(vcfvar)
            vcf_file.close()
        except IOError:
            self.vaselogger.warning(f"Could not open VCF file {vcf_fileloc}")
        finally:
            return sample_variant_list

    def determine_variant_type(self, vcfvariantstart, vcfvariantstop):
        """Determines and returns the variant type.

        Determination of the variant type is based on the distance between the right and left most genomic position of
        the variant reference and alternative allele(s).

        Parameters
        ----------
        vcfvariantstart : int
            Leftmost genomic position of the variant
        vcfvariantstop : int
            Rightmost genomic position of the variant

        Returns
        -------
        str
            Type of variant (snp/indel)
        """
        if (vcfvariantstop - vcfvariantstart) <= 1:
            return "snp"
        return "indel"

    def determine_read_search_window(self, varianttype, vcfvariant):
        """Determines and returns the search window for fetching reads.

        The determined search window depends on the provided variant type. SNPs return a search window of SNP position
        -1 and +1. Indels return a search window based indel lefmost position and length. If the variant type is
        neither SNP or indel, the search window returned will be [-1, -1].

        Parameters
        ----------
        varianttype : str
            Variant type to determine search window for
        vcfvariant : VcfVariant
            Variant to determine search window with

        Returns
        -------
        list of int
            Search window start and stop, -1 and -1 if variant type is invalid
        """
        if varianttype == "snp":
            return [vcfvariant.get_variant_pos() - 1, vcfvariant.get_variant_pos() + 1]
        elif varianttype == "indel":
            return self.determine_indel_read_range(vcfvariant.get_variant_pos(),
                                                   vcfvariant.get_variant_ref_allele(),
                                                   vcfvariant.get_variant_alt_alleles())
        return [-1, -1]

    # Returns the search start and stop to use for searching BAM reads overlapping with the range of the indel.
    def determine_indel_read_range(self, variantpos, variantref, variantalts):
        """Determines and returns the search start and stop to use for an indel.

        Parameters
        ----------
        variantpos : int
            Leftmost genomic position of the variant
        variantref : str
            Variant reference allele(s)
        variantalts : tuple
            Variant alternative allele(s)

        Returns
        -------
        list of int
            Search window start and stop
        """
        searchstart = variantpos
        searchstop = variantpos + max(
            max([len(x) for x in variantref.split(",")]),
            max([len(x) for x in variantalts])
        )
        return [searchstart, searchstop]

    def get_variant_reads(self, contextid, variantchrom,
                          variantstart, variantend,
                          bamfile,
                          write_unm=False, umatelist=None):
        """Fetches and returns reads overlapping with a specified variant.

        First reads overlapping directly with the variant position are fetched. Then the read mates are fetched using
        the RNEXT and PNEXT values of each read. Lastly, it is ensured that each read only occurs once.

        Parameters
        ----------
        contextid : str
            Identifier of the context to fetch reads for
        variantchrom : str
            Chromosome name to fetch reads from
        variantstart : int
            Leftmost genomic position to use for fetching
        variantend : int
            Rightmost genomic position to use for fetching
        bamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile
        write_unm : bool
            Save read identifiers with an unmapped mate
        umatelist : list of str
            Identifiers of reads with an unmapped mate

        Returns
        -------
        variantreads : list of pysam.AlignedSegment
            Fetched reads and their read mates
        """
        hardclipped_read_num = 0
        duplicate_read_num = 0
        secondary_read_num = 0

        read_objects = []
        variantreads = []
        clipped_reads = []
        list_r1 = []
        list_r2 = []

        for vread in bamfile.fetch(variantchrom, variantstart, variantend):
            if vread.is_duplicate:
                duplicate_read_num += 1
            if vread.is_secondary:
                secondary_read_num += 1
            if self.read_is_hard_clipped(vread):
                hardclipped_read_num += 1
                clipped_reads.append(vread)
                continue
            read_objects.append(vread)

        self.vaselogger.debug(f"Fetched {hardclipped_read_num} reads with hardclipped bases")
        for clipped in clipped_reads:
            read_objects.append(self.fetch_primary_from_secondary(clipped, bamfile))

        for vread in read_objects:
            if vread.is_read1:
                list_r1.append(vread)
            elif vread.is_read2:
                list_r2.append(vread)

        list_r1_ids = [x.query_name for x in list_r1]
        list_r2_ids = [x.query_name for x in list_r2]

        for r1 in list_r1:
            if r1.query_name not in list_r2_ids:
                list_r2.append(self.fetch_mate_read(r1.query_name, r1.next_reference_name, r1.next_reference_start,
                                                    self.get_read_pair_num(r1), bamfile))
        for r2 in list_r2:
            if r2.query_name not in list_r1_ids:
                list_r1.append(self.fetch_mate_read(r2.query_name, r2.next_reference_name, r2.next_reference_start,
                                                    self.get_read_pair_num(r2), bamfile))

        variantreads = list_r1 + list_r2
        variantreads = self.uniqify_variant_reads(variantreads)
        return variantreads

    def fetch_primary_from_secondary(self, secondary_read, bamfile):
        """Fetches and returns the primary alignment of a read based on the position recorded in its SA tag.

        The SA tag is read from a provided Pysam read object. Then, reads at the recorded alignment position are
        fetched. The first non-secondary alignment with a matching read ID and pair number is returned.

        Parameters
        ----------
        secondary_read : pysam.AlignedSegment
            Secondary alignment whose primary alignment will be located
        bamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile

        Returns
        -------
        primary_read : pysam.AlignedSegment
            Primary alignment with the same read ID and pair number as the input secondary alignment.
        """
        primary_locus = secondary_read.get_tag("SA").split(",")[0:2]
        for fetched in bamfile.fetch(primary_locus[0], int(primary_locus[1]) - 1, int(primary_locus[1])):
            if ((fetched.query_name == secondary_read.query_name)
                    and (fetched.is_read1 == secondary_read.is_read1)
                    and not fetched.is_secondary):
                primary_read = fetched
                break
        return primary_read

    def uniqify_variant_reads(self, variantreads):
        """Ensures each read only occurs once and returns the updates set.

        Parameters
        ----------
        variantreads : list of pysam.AlignedSegment
            Sets of reads to process

        Returns
        -------
        unique_variantreads : list of pysam.AlignedSegment
            Set of reads with each read occuring only once
        """
        unique_variantreads = []
        checklist = []
        for fetched in variantreads:
            readpn = "2"
            if fetched.is_read1:
                readpn = "1"
            id_pair = (fetched.query_name, readpn)
            if id_pair not in checklist:
                unique_variantreads.append(fetched)
                checklist.append(id_pair)
        return unique_variantreads

    def fetch_mate_read(self, readid, rnext, pnext, pair_num, bamfile):
        """Fetches and returns the mate read of a specified read.

        The mate read is fetched from an opened alignment file using the RNEXT and PNEXT value of the read. Of the
        fetched reads, only the read with the same identifier is returned.

        Parameters
        ----------
        rnext : str
            Chromosome name the mate read is located on
        pnext : int
            Leftmost genomic position of the mate read
        readid : str
            Identifier of the read to search the mate for
        bamfile : pysam.AlignmentFile
            Already openend pysam alignment file

        Returns
        -------
        pysam.AlignedSegment or None
            The mate read if found, None if not
        """
        for bamread in bamfile.fetch(rnext, pnext, pnext + 1):
            if (bamread.query_name == readid
                    and bamread.reference_start == pnext
                    and self.get_read_pair_num(bamread) != pair_num):
                return bamread
        return None

    def determine_context(self, contextreads, contextorigin, contextchr):
        """Determines and returns an acceptor/donor context.

        The acceptor/donor context is determined from a set of reads overlapping with the variant including their read
        mates. The leftmost and rightmost genomic positions of all reads are collected and filtered for outliers.
        The context start (leftmost genomic) position is then determined by taking minimum and maximum leftmost and
        rightmost position respectively.

        Parameters
        ----------
        contextreads: list of pysam.AlignedSegment
        contextorigin : int
            Variant genomic position the context will be based on
        contextchr : str
            Chromosomae name the context is located on

        Returns
        -------
        list of str and int
            Essential context data (chromosome, variant pos, start, end)
        """
        # Check whether there are reads to determine the context for.
        if not contextreads:
            return []

        # Get read start and stop position, while filtering out read mates that map to different chr.
        starts = [conread.reference_start
                  for conread in contextreads
                  if conread.reference_name == contextchr]

        stops = [conread.reference_end
                 for conread in contextreads
                 if conread.reference_name == contextchr and conread.reference_end is not None]

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

        self.vaselogger.debug(f"{num_diff_chr} read mate(s) filtered due to alignment to different reference sequence.")
        self.vaselogger.debug(f"{num_filtered} outlier read position(s) filtered.")

        # Set variant context as chr, min start, max end.
        contextstart = min(filtered_starts)
        contextend = max(filtered_stops)
        return [contextchr, contextorigin, contextstart, contextend]

    def filter_outliers(self, pos_list, k=3):
        """Filters outliers from a list of start or stop positions and returns the filtered list.

        Outlier start/stop positions are filtered from the list using Tukey's Fences method. For more info please see
        https://en.wikipedia.org/wiki/Outlier#Tukey's_fences

        Parameters
        ----------
        pos_list : list of int
            Start/stop positions
        k : int
            Factor to determine outlier
        Returns
        -------
        """
        # First and third quartile values of the positions.
        q1 = np.percentile(pos_list, 25)
        q3 = np.percentile(pos_list, 75)
        # Interquartile range.
        iq = q3 - q1
        # Only include positions that fall within the range of (q1 to q3) +/- k*iq.
        filtered = [x for x in pos_list
                    if ((q1 - (k * iq)) <= x <= (q3 + (k * iq)))]
        return filtered

    def determine_largest_context(self, contextorigin, acceptor_context,
                                  donor_context):
        """Determines the size of the variant context based on both the acceptor and donor reads.

        Parameters
        ----------
        contextorigin : int
            Variant position that the context is based on
        acceptor_context : list of str and int
            Essential acceptor context data
        donor_context : list of str and int
            Essential donor context data

        Returns
        -------
        list of str and int
            Context (chromosome, variant pos, start, end)
        """
        largest_context = [donor_context[0]]
        largest_context.append(contextorigin)
        # Determine and save the smallest context start.
        largest_context.append(min(acceptor_context[2], donor_context[2]))
        # Determine and save the largest context end.
        largest_context.append(max(acceptor_context[3], donor_context[3]))
        return largest_context

    # ===METHODS TO PRODUCE THE VALIDATION FASTQ FILES=========================
    def build_donor_fq(self, donorbamreaddata, fr, fastq_outpath):
        """Build and writes a fastq file containing only donor reads.

        Parameters
        ----------
        donorbamreaddata : list of tuple
            Donor reads to write to fastq file
        fr : str
            Write forward ('1') or reverse ('2') fastq file
        fastq_outpath : str
            Path and name to write donor fastq file to
        """
        # Write the new VaSe FastQ file.
        # vasefq_outname = self.set_fastq_out_path(fastq_outpath, fr, 1)
        fqgz_outfile = io.BufferedWriter(open(fastq_outpath, "wb"))

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

    def is_required_read(self, bamread, fr):
        """Checks and returns whether the current read is the one that is required.

        When writing the validation fastq files only R1, or forward (F), reads should go to the R1 validation fastq
        file. Similarly, the R2 reads should only go into the R2 validation fastq file.

        Parameters
        ----------
        bamread : DonorBamRead
            Read to check
        fr : str
            Whether the read is required for R1 or R2

        Returns
        -------
        bool
            True if the read is required, False if not
        """
        if fr == "F":
            return bamread.is_read1
        return bamread.is_read2

    def set_fastq_out_path(self, outpath, fr, lnum):
        """Sets and returns the fastq output path and filename.

        Parameters
        ----------
        outpath : str
            Path and suffix to write fastq file to
        fr: str
            Will the fastq file be R1 ('1') or R2 ('2')
        lnum: int
            Lane number (i.e.: 1, 2, 3 or 4)

        Returns
        -------
        str
            Full path to write fastq file to
        """
        if fr == "1":
            return f"{outpath}_{datetime.now().date()}_L{lnum}_R1.fastq"
        return f"{outpath}_{datetime.now().date()}_L{lnum}_R2.fastq"

    # ===METHODS TO OBTAIN SOME DATA OF THE VASEBUILDER OBJECT=================
    def get_creation_id(self):
        """Returns the identifier of the current VaSeBuilder object.

        Returns
        -------
        self.creation_id : str
            VaSeBuilder creation identifier
        """
        return self.creation_id

    def get_creation_time(self):
        """Returns the date and time the current VaSeBuilder object has been made.

        Returns
        -------
        str
            VaSeBuilder creation date and time
        """
        return self.creation_time

    def get_vcf_variant_id(self, vcfvariant):
        """Returns an identifier for a variant as CHROM_POS.

        Parameters
        ----------
        vcfvariant : pysam.VariantRecord
            Variant for which to construct an identifier
        Returns
        -------
        str
            Variant identfier as CHROM_POS
        """
        return f"{vcfvariant.chrom}_{vcfvariant.pos}"

    def get_read_pair_num(self, pysam_bamread):
        """Returns the read pair number of a provided read.

        Parameters
        ----------
        pysam_bamread : pysam.AlignedSegment
            Read to determine the read pair number of

        Returns
        -------
        str
            Read pair number ('1' or '2')
        """
        if pysam_bamread.is_read1:
            return "1"
        return "2"

    # ===METHODS TO WRITE OUTPUT FILES=========================================
    def write_used_donor_files(self, outfileloc, filesamplemap,
                               used_donor_files, vbuuid):
        """Writes the donor alignment or variant files used in constructing variant contexts to an output file.

        Parameters
        ----------
        outfileloc : str
            Path to write the output file to
        filesamplemap : dict
            Donor files per sample
        used_donor_files : list
            Donor alignment or variant files used to construct varian contexts
        vbuuid : str
            Unique identifier of the current VaSeBuilder
        """
        try:
            with open(outfileloc, "w") as outfile:
                outfile.write(f"#VBUUID: {vbuuid}\n")
                outfile.write("#SampleId\tDonorFile\n")
                for sampleid, samplefile in filesamplemap.items():
                    if samplefile in used_donor_files:
                        outfile.write(f"{sampleid}\t{samplefile}\n")
        except IOError as ioe:
            self.vaselogger.critical("Could not write used donor files to "
                                     f"{outfileloc}")

    def write_optional_output_files(self, outpath, contextfile):
        """Writes optional output files to a specified output location.

        The optional output files contain data about acceptor and donor contexts as well as output files containing the
        reads with unmapped mates per variant context. The optional output files are only produced when VaSeBuilder is
        run with with debug set to True.

        Parameters
        ----------
        outpath: str
            Path to write the optional output file to
        contextfile: VariantContextFile
            Variant context data
        """
        # Write the optional acceptor context files; acceptor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaselogger.debug("Writing acceptor contexts to "
                              f"{outpath}acceptorcontexts.txt")
        contextfile.write_acceptor_context_file(f"{outpath}acceptorcontexts.txt", self.creation_id)

        self.vaselogger.debug("Writing acceptor context statistics to "
                              f"{outpath}acceptorcontextstats.txt")
        contextfile.write_acceptor_context_stats(f"{outpath}acceptorcontextstats.txt", self.creation_id)

        self.vaselogger.debug("Writing acceptor context read identifiers with "
                              "unmapped mates to"
                              f"mates to {outpath}acceptor_unmapped.txt")
        contextfile.write_acceptor_unmapped_mates(f"{outpath}acceptor_unmapped.txt", self.creation_id)

        self.vaselogger.debug("Writing left and right most read positions of "
                              "each acceptor context to "
                              f"{outpath}acceptor_positions.txt")
        contextfile.write_acceptor_left_right_positions(f"{outpath}acceptor_positions.txt", self.creation_id)

        # Write the optional donor context files; donor contexts, read ids with unmapped mate and left/right positions.
        self.vaselogger.debug("Writing donor contexts to "
                              f"{outpath}donorcontexts.txt")
        contextfile.write_donor_context_file(f"{outpath}donorcontexts.txt", self.creation_id)

        self.vaselogger.debug("Writing donor context statistics to "
                              f"{outpath}donorcontextstats.txt")
        contextfile.write_donor_context_stats(f"{outpath}donorcontextstats.txt", self.creation_id)

        self.vaselogger.debug("Writing donor context read identifiers "
                              f"with unmapped mates to {outpath}donor_unmapped.txt")
        contextfile.write_donor_unmapped_mates(f"{outpath}donor_unmapped.txt", self.creation_id)

        self.vaselogger.debug("Writing left and right most read positions "
                              f"of each donor context to {outpath}donor_positions.txt")
        contextfile.write_donor_left_right_positions(f"{outpath}donor_positions.txt", self.creation_id)

        # Write the optional variant context files; acceptor & donor unmapped mates and left/right positions.
        self.vaselogger.debug("Writing variant context acceptor read "
                              f"identifiers with unmapped mates to {outpath}varcon_unmapped_acceptor.txt")
        contextfile.write_reads_with_unmapped_mate("acceptor", f"{outpath}varcon_unmapped_acceptor.txt",
                                                   self.creation_id)

        self.vaselogger.debug("Writing variant context donor read identifiers "
                              f"with unmapped mates to {outpath}varcon_unmapped_donor.txt")
        contextfile.write_reads_with_unmapped_mate("donor", f"{outpath}varcon_unmapped_donor.txt", self.creation_id)

        self.vaselogger.debug("Writing variant context left and right most "
                              f"read positions of acceptor reads to {outpath}varcon_positions_acceptor.txt")
        contextfile.write_left_right_positions(
            "acceptor",
            f"{outpath}varcon_positions_acceptor.txt", self.creation_id
        )

        self.vaselogger.debug("Writing variant context left and right most "
                              "read positions of donor reads to "
                              f"{outpath}varcon_positions_donor.txt")
        contextfile.write_left_right_positions(
            "donor",
            f"{outpath}varcon_positions_donor.txt", self.creation_id
        )

    def write_bed_file(self, variantcontextdata, bedoutloc, vbuuid):
        """Writes variant contexts as a BED file

        Parameters
        ----------
        variantcontextdata : list of VariantContext
            Variants contexts to write to BED file
        bedoutloc : str
            Path to write the BED file to
        """
        try:
            with open(bedoutloc, "w") as bedoutfile:
                bedoutfile.write(f"#VBUUID: {vbuuid}\n")
                for varcon in variantcontextdata:
                    bedoutfile.write(f"{varcon.get_variant_context_chrom()}\t{varcon.get_variant_context_start()}\t"
                                     f"{varcon.get_variant_context_end()}\t{varcon.get_variant_context_id()}\n")
        except IOError:
            self.vaselogger.warning(f"Could not write variant context data to BED file: {bedoutloc}")

    def check_sequence_names(self, referencefile, alignmentfile):
        """Checks and returns whether the chromosome names in the genome reference and alignment file are the same.

        Parameters
        ----------
        referencefile:
        alignmentfile:

        Returns
        -------

        """
        reference_seqnames = self.get_reference_sequence_names(referencefile)
        alignment_seqnames = self.vb_scanner.get_alignment_sequence_names(alignmentfile)
        shared_seqnames = reference_seqnames & alignment_seqnames
        if len(shared_seqnames) < len(reference_seqnames) or len(shared_seqnames) < len(alignment_seqnames):
            self.vaselogger.warning("Reference file and alignment file do not contain the same sequence names")

    def get_reference_sequence_names(self, reference_fileloc):
        """Returns the sequence names from the reference genome fasta file

        Parameters
        ----------
        reference_fileloc: str
            Path to genomic reference fasta file

        Returns
        -------
        reference_seqnames : list of str
            Extracted chromosome names
        """
        reference_seqnames = set()
        try:
            with open(reference_fileloc, "r") as reference_file:
                for fileline in reference_file:
                    if fileline.startswith(">"):
                        filelinedata = fileline.strip().split(" ")
                        reference_seqnames.add(filelinedata[0][1:])
        except IOError:
            self.vaselogger.critical(f"Could not read genome reference file {reference_fileloc}")
            exit()
        return reference_seqnames

    def divide_donorfastqs_over_acceptors(self, list_of_donorfqs, num_of_acceptors):
        """Divides a list of donor fastq files over a number of acceptor fastq files.

        Donor fastq files are divided equally into a specified number of groups. The number of groups is specified by
        the number of acceptor fastq files to use as template.

        Parameters
        ----------
        list_of_donorfqs: list of str
            Paths to donor fastq files
        num_of_acceptors: int
            Number of template fastq files
        """
        cycler = 0
        split_donors = [[] for x in range(num_of_acceptors)]
        copy_list_donorfqs = list(list_of_donorfqs)
        while copy_list_donorfqs:
            cycle = cycler % num_of_acceptors
            split_donors[cycle].append(copy_list_donorfqs.pop())
            cycler += 1
        return split_donors

    # BUILDS A SET OF R1/R2 VALIDATION FASTQS WITH ALREADY EXISTING DONOR FASTQS
    def run_d_mode(self, variantcontextfile, fq_out):
        """Runs VaSeBuilder D-mode.

        This run mode produces one R1 and R2 file with donor reads of all variant contexts. This mode does not produce
        validation fastq files.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Variants context to use
        fq_out : str
            Path and suffix to write donor fastq files to
        """
        self.vaselogger.info("Running VaSeBuilder D-mode")
        # Combine all donor reads from all variant contexts
        add_list = variantcontextfile.get_all_variant_context_donor_reads()
        self.vaselogger.info("Writing FastQ files.")
        # Write the donor fastq files.
        for i in ["1", "2"]:
            self.vaselogger.info(f"Start writing the R{i} FastQ files.")
            fq_starttime = time.time()
            fqoutpath = self.set_fastq_out_path(fq_out, i, 1)
            self.build_donor_fq(add_list, i, fqoutpath)
        self.vaselogger.info(f"Wrote all R{i} FastQ files.")
        self.vaselogger.debug(f"Writing R{i} FastQ file(s) took {time.time() - fq_starttime} seconds.")
        self.vaselogger.info("Finished writing donor FastQ files.")

    def run_f_mode(self, variantcontextfile, fq1_in, fq2_in, fq_out, random_seed):
        """Runs VaSeBuilder F-mode.

        This run mode creates a full set of validation fastq files from established variant contexts and a set of
        acceptor fastq files to use as template to create the validation fastq files.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Variant contexts to use
        fq1_in : list of str
            R1 fastq files to use as template
        fq2_in : list of str
            R2 fastq files to use as template
        fq_out : str
            Path and suffix to write validation fastq files to
        """
        self.vaselogger.info("Running VaSeBuilder F-mode")

        # Combine all donor reads from all variant contexts
        add_list = variantcontextfile.get_all_variant_context_donor_reads()
        self.vaselogger.info("Writing FastQ files.")
        skip_list = set(variantcontextfile.get_all_variant_context_acceptor_read_ids())
        # skip_list = list(set(variantcontextfile.get_all_variant_context_acceptor_read_ids()))

        donor_read_add_data = {}
        for i, fq_i in zip(["1", "2"], [fq1_in, fq2_in]):
            # Write the fastq files.
            self.vaselogger.info(f"Start writing the R{i} FastQ files.")
            fq_starttime = time.time()
            self.build_fastq_v2(fq_i, skip_list, add_list, i, fq_out, random_seed, donor_read_add_data)
            self.vaselogger.info(f"Wrote all R{i} FastQ files.")
            self.vaselogger.debug(f"Writing R{i} FastQ file(s) took {time.time() - fq_starttime} seconds.")
        self.vaselogger.info("Finished writing FastQ files.")
        self.write_donor_insert_positions_v2(donor_read_add_data, f"{fq_out}_donor_read_insert_positions.txt")

    def run_p_mode(self, variantcontextfile, outpath, fq_out):
        """Runs VaSeBuilder P-mode.

        This run mode produces an R1 and R2 fastq file with donor reads for each variant context. This mode does not
        create or write validation fastq files.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Established variant contexts
        outpath : str
            Path to folder to write output files to
        fq_out : str
            Path and suffix to write fastq out files to
        """
        context_fq_links = {}

        self.vaselogger.info("Running VaSeBuilder P-mode")
        self.vaselogger.info("Begin writing variant FastQ files.")
        variantcontexts = variantcontextfile.get_variant_contexts()
        for context in variantcontexts:
            add_list = context.get_donor_read_strings()
            self.vaselogger.debug(f"Writing variant FastQs for variant {context.context_id}.")

            r1_donorfq = self.set_fastq_out_path(fq_out + context.context_id, "1", 1)
            r2_donorfq = self.set_fastq_out_path(fq_out + context.context_id, "2", 1)
            self.build_donor_fq(add_list, "1", r1_donorfq)
            self.build_donor_fq(add_list, "2", r2_donorfq)
            context_fq_links[context.context_id] = [r1_donorfq, r2_donorfq]
        self.vaselogger.info("Finished writing variant FastQ files.")

        # Write the P-mode link file (links context to fastq files for the current run)
        self.write_pmode_linkfile(outpath, context_fq_links)

    def run_x_mode(self, sampleidlist, donorvcfs, donorbams, acceptorbam, genomereference, outdir, varconout,
                   variantlist):
        """Runs VaSeBuilder X-mode.

        This run mode only produces the variant contexts and writes the associated output files. This mode does not
        create or write validation fastq files.

        Parameters
        ----------
        sampleidlist : list of str
            Names/identifiers of samples to process
        donorvcfs : dict
            Donor variant files per sample
        donorbams : dict
            Donor alignment files per sample
        acceptorbam : str
            Path to alignment file to use as acceptor
        genomereference : str
            Path to genome reference fasta file
        outdir : str
            Path to write output files to
        varconout : str
            Path to folder to write the variant context file to
        variantlist : dict
            Variants to process per sample
        """
        self.vaselogger.info("Running VaSeBuilder X-mode")
        self.bvcs(sampleidlist, donorvcfs, donorbams, acceptorbam, outdir, genomereference, varconout, variantlist)

    # =====SPLITTING THE BUILD_VARCON_SET() INTO MULTIPLE SMALLER METHODS=====
    def bvcs(self, sampleidlist, vcfsamplemap, bamsamplemap, acceptorbamloc, outpath, reference_loc, varcon_outpath,
             variantlist):
        """Builds and returns a variant context file.

        The variant context file is build from the provided variants and acceptor and donors alignment files.

        Parameters
        ----------
        sampleidlist : list of str
            Sample name/identifier
        vcfsamplemap : dict
            Variant files per sample
        bamsamplemap : dict
            Alignment files per sample
        acceptorbamloc : str
            Path to alignment file to use as acceptor
        outpath : str
            Path to folder to write output files to
        reference_loc : str
            Path to the genomic reference fasta file
        varcon_outpath : str
            Path and name to write variant context file to
        variantlist : dict
            Variants to use per sample

        Returns
        -------
        variantcontexts : VariantContextFile
            Established variant contexts
        """
        donor_vcfs_used = []
        donor_bams_used = []
        variantcontexts = VariantContextFile()

        try:
            acceptorbamfile = pysam.AlignmentFile(acceptorbamloc, reference_filename=reference_loc)
        except IOError:
            self.vaselogger.critical("Could not open Acceptor BAM/CRAM")
            exit()

        # Start iterating over the samples
        for sampleid in sampleidlist:
            self.vaselogger.debug(f"Start processing sample {sampleid}")
            sample_variant_filter = self.bvcs_set_variant_filter(sampleid, variantlist)
            samplevariants = self.get_sample_vcf_variants(vcfsamplemap[sampleid], sample_variant_filter)

            if not samplevariants:
                self.vaselogger.warning(f"No variants obtained for sample {sampleid}. Skipping sample")
                continue

            # Call the method that will process the sample
            self.bvcs_process_sample(sampleid, variantcontexts, acceptorbamfile, bamsamplemap[sampleid],
                                     reference_loc, samplevariants)

            # Add the used donor VCF and BAM to the lists of used VCF and BAM files
            donor_bams_used.append(vcfsamplemap[sampleid])
            donor_vcfs_used.append(bamsamplemap[sampleid])

        # Check if there are no variant contexts.
        if variantcontexts.get_number_of_contexts() <= 0:
            self.vaselogger.info("No variant contexts were created. No output files will be written.")
            return None

        # Write the output variant context output data
        self.bvcs_write_output_files(outpath, varcon_outpath, variantcontexts, vcfsamplemap, bamsamplemap,
                                     donor_vcfs_used, donor_bams_used)

        # Checks whether the program is running on debug and, if so, write some extra output files.
        if self.vaselogger.getEffectiveLevel() == 10:
            self.write_optional_output_files(outpath, variantcontexts)
        return variantcontexts

    def bvcs_process_sample(self, sampleid, variantcontextfile, abamfile, dbamfileloc, referenceloc, samplevariants):
        """Processes a sample and adds variant contexts to a variant context file.

        Parameters
        ----------
        sampleid : str
            Sample name/identifier
        variantcontextfile : VariantContextFile
            Variant context file that saves variant contexts
        abamfile: pysam.AlignmentFile
            Already opened pysam AlignmentFile to use as acceptor
        dbamfileloc: str
            Path to the alignment file to use as donor
        referenceloc: str
            Path to the genomic reference fasta file
        samplevariants : list of VcfVariants
            Variants to process for the specified sample
        """
        try:
            donorbamfile = pysam.AlignmentFile(dbamfileloc, reference_filename=referenceloc)
        except IOError:
            self.vaselogger.warning(f"Could not open {dbamfileloc} ; Skipping {sampleid}")
            return

        # Iterate over the sample variants
        for samplevariant in samplevariants:
            variantcontext = self.bvcs_process_variant(sampleid, variantcontextfile, samplevariant, abamfile,
                                                       donorbamfile, True)
            if not variantcontext:
                self.vaselogger.info(f"Could not establish variant context ; Skipping.")
                continue

            varcon_collided = variantcontextfile.context_collision_v2(variantcontext.get_context())
            if varcon_collided is not None:
                self.vaselogger.debug(
                    f"Variant context {variantcontext.get_variant_context_id()} overlaps with variant context"
                    f"{varcon_collided.get_variant_context_id()}")
                if varcon_collided.get_variant_context_sample() != variantcontext.get_variant_context_sample():
                    self.vaselogger.debug(f"Colliding contexts from different sample ; Skipping.")
                    continue
                self.vaselogger.debug(f"Merging variant contexts {variantcontext.get_variant_context_id()} with "
                                      f"{varcon_collided.get_variant_context_id()}")
                variantcontext = self.merge_variant_contexts(varcon_collided, variantcontext)
            variantcontextfile.add_existing_variant_context(variantcontext.get_variant_context_id(), variantcontext)

            # If this widest context overlaps an existing variant context, skip it.
            #if variantcontextfile.context_collision(variantcontext.get_context()):
            #    self.vaselogger.debug(f"Variant context {variantcontext.get_variant_context_id()} overlaps with an"
            #                          f"already existing variant context; Skipping.")
            #    continue
            #variantcontextfile.add_existing_variant_context(variantcontext.get_variant_context_id(), variantcontext)

    def bvcs_process_variant(self, sampleid, variantcontextfile, samplevariant, abamfile, dbamfile, write_unm=False):
        """Processes a variant and returns the established variant context.

        Parameters
        ----------
        sampleid : str
            Sample name/identifier
        variantcontextfile: VariantContextFile
            Variant context file that saves the variant contexts
        samplevariant : pysam.VariantRecord
            Variant to process and for which to establish an accpetor, donor and variant context
        abamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile to use as acceptor
        dbamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile to use as donor
        write_unm : bool
            Whether to save read identifiers with unmapped read mates

        Returns
        -------
        vcontext : VariantContext or None
            The established variant context, None if context could not be established
        """
        # Establish the donor context
        variantid = self.get_vcf_variant_id(samplevariant)
        varianttype = self.determine_variant_type(samplevariant.start, samplevariant.stop)

        self.vaselogger.debug(f"Processing variant {variantid}.")
        self.debug_msg("vw", variantid)
        self.vaselogger.debug(f"Variant {variantid} determined to be {varianttype}")
        searchwindow = [samplevariant.start, samplevariant.stop]
        self.vaselogger.debug(f"Search window determined to be {samplevariant.chrom}:"
                              f"{searchwindow[0]}-{searchwindow[1]}")

        # Check whether the variant overlaps with an existing variant context
        overlapping_varcon = variantcontextfile.variant_is_in_context(varianttype, samplevariant.chrom, *searchwindow)
        if overlapping_varcon is not None:
            self.vaselogger.debug(f"Variant {variantid} is located in an existing context")
            if sampleid != overlapping_varcon.get_variant_context_sample():
                self.vaselogger.debug(f"Overlapping variant and variant context are not of the same sample ; "
                                      f"Skipping variant.")
                return None

        # Determine the donor context.
        self.debug_msg("dc", variantid)
        t0 = time.time()
        dcontext = self.bvcs_establish_context(sampleid, variantid, samplevariant.chrom, samplevariant.start,
                                               searchwindow, dbamfile, write_unm)
        if not dcontext:
            self.vaselogger.info(f"Could not establish donor context ; Skipping variant {variantid}")
            return None
        self.debug_msg("dc", variantid, t0)
        self.vaselogger.debug(f"Donor context determined to be {dcontext.get_context_chrom()}:"
                              f"{dcontext.get_context_start()}-{dcontext.get_context_end()}")

        # Determine the acceptor context.
        self.debug_msg("ac", variantid)
        t0 = time.time()
        acontext = self.bvcs_establish_context(sampleid, variantid, samplevariant.chrom, samplevariant.start,
                                               searchwindow, abamfile, write_unm, dcontext.get_context())
        self.debug_msg("ac", variantid, t0)
        self.vaselogger.debug(f"Acceptor context determined to be {acontext.get_context_chrom()}:"
                              f"{acontext.get_context_start()}-{acontext.get_context_end()}")

        # Determin the variant context.
        self.debug_msg("cc", variantid)
        t0 = time.time()
        vcontext = self.bvcs_establish_variant_context(sampleid, variantcontextfile, variantid, samplevariant.chrom,
                                                       samplevariant.start, acontext, dcontext, abamfile, dbamfile,
                                                       write_unm)
        self.debug_msg("cc", variantid, t0)
        if vcontext is not None:
            self.vaselogger.debug(f"Combined context determined to be {vcontext.get_variant_context_chrom()}:"
                                  f"{vcontext.get_variant_context_start()}-{vcontext.get_variant_context_end()}")
        return vcontext

    def bvcs_establish_variant_context(self, sampleid, variantcontextfile, variantid, variantchrom, variantpos,
                                       acontext, dcontext, abamfile, dbamfile, write_unm=False):
        """Establishes and returns a variant context.

        The variant context window is established based on the acceptor and donor contexts. Reads and mates overlapping
        with the variant context are also fetched and added to the variant context.

        Parameters
        ----------
        sampleid : str
            Sample name/identifier
        variantcontextfile : VariantContextFile
            Variant context file that saves variant contexts
        variantid : str
            Variant identifier
        variantchrom : str
            Chromosome the
        variantpos: int
            Leftmost genomic position of the variant
        acontext : OverlapContext
            Established acceptor context
        dcontext : OverlapContext
            Established donor context
        abamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile used as acceptor
        dbamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile used as donor
        write_unm : bool
            Whether to save identifiers of reads with unmapped mates

        Returns
        -------
        variant_context : VariantContext
            Established variant context with all associated data
        """
        unmapped_alist = []
        unmapped_dlist = []

        # Determine the variant context from the acceptor and donor context
        vcontext_window = self.determine_largest_context(variantpos, acontext.get_context(), dcontext.get_context())
        #if variantcontextfile.context_collision(vcontext_window):
        #    self.vaselogger.debug(f"Variant context {variantid} overlaps with an already existing context; Skipping.")
        #    return None

        # Gather variant context donor reads.
        self.debug_msg("cdr", variantid)
        t0 = time.time()
        vcontext_dreads = self.get_variant_reads(variantid, vcontext_window[0], vcontext_window[2], vcontext_window[3],
                                                 dbamfile, write_unm, unmapped_dlist)
        self.debug_msg("cdr", variantid, t0)

        # Gather variant context acceptor reads.
        self.debug_msg("car", variantid)
        t0 = time.time()
        vcontext_areads = self.get_variant_reads(variantid, vcontext_window[0], vcontext_window[2], vcontext_window[3],
                                                 abamfile, write_unm, unmapped_alist)
        self.debug_msg("car", variantid, t0)

        variant_context = VariantContext(variantid, sampleid, *vcontext_window, vcontext_areads, vcontext_dreads,
                                         acontext, dcontext)

        # Set the variant context unmapped read ids
        variant_context.set_unmapped_acceptor_mate_ids(unmapped_alist)
        variant_context.set_unmapped_donor_mate_ids(unmapped_dlist)
        return variant_context

    # Establishes an acceptor/donor context by fetching reads (and their mates) overlapping directly with the variant
    def bvcs_establish_context(self, sampleid, variantid, variantchrom, variantpos, searchwindow, bamfile,
                               write_unm=False,
                               fallback_window=None):
        """Establishes and returns an acceptor/donor context.

        The context window is established from fetched reads and their mates overlapping with the variant.

        Parameters
        ----------
        sampleid : str
            Sample name/identifier
        variantid : str
            Variant indentifier
        variantchrom : str
            Chromosome name the variant is located on
        variantpos : int
            Leftmost genomic position of the variant
        searchwindow: list of int
            Sreach window to use for fetching reads
        bamfile : pysam.AlignmentFile
            Already opened pysam AlignmentFile
        write_unm: bool
            Save the read identifiers with unmapped mates
        fallback_window : list
            Window to set if no window could be set/determined

        Returns
        -------
        OverlapContext or None
            Acceptor/donor context if context is established, None if not
        """
        unmappedlist = []

        # Fetch context reads and establish the context window
        self.vaselogger.debug("Fetching reads.")
        t0 = time.time()
        context_reads = self.get_variant_reads(variantid, variantchrom, *searchwindow, bamfile, write_unm, unmappedlist)
        self.vaselogger.debug(f"Fetching reads took {time.time() - t0} seconds")

        if not context_reads:
            self.vaselogger.debug("No reads were found.")
        context_window = self.determine_context(context_reads, variantpos, variantchrom)

        # Check whether the context_window is valid and, if not, whether a fallback window has been set
        if not context_window:
            self.vaselogger.debug("Context window could not be determined.")
            if fallback_window:
                context_window = fallback_window
                self.vaselogger.debug("Context window is set to the fall back window.")
            else:
                return None

        # Create the context object containing all context data
        adcontext = OverlapContext(variantid, sampleid, *context_window, context_reads)
        adcontext.set_unmapped_mate_ids(unmappedlist)
        return adcontext

    # Sets the variant list for a sample if applicable, otherwise return None
    def bvcs_set_variant_filter(self, sampleid, variant_list):
        """Sets the variant list for a sample if applicable, otherwise return None

        Parameters
        ----------
        sampleid : str
            Sample name/identifier
        variant_list : dict
            Variant tuples with chromosome name and genomic position per sample

        Returns
        -------
        list of tuple
            Variants tuples with chromosome name and position
        """
        sample_variant_filter = None
        if variant_list is not None:
            if sampleid in variant_list:
                sample_variant_filter = variant_list[sampleid]
            else:
                sample_variant_filter = []
        return sample_variant_filter

    def bvcs_write_output_files(self, outpath, varcon_outpath, variantcontextfile, vcfsamplemap, bamsamplemap,
                                used_donor_vcfs, used_donor_bams):
        """Writes VaSeBuilder output files.

        These output files are the standard output files when building a variant context file. Output files include
        the variant context file, variant context statistics file, variant context BED file, used donor variant files
        and used donor alignment files.

        Parameters
        ----------
        outpath: str
            Path to folder to write the output files to
        varcon_outpath : str
            Path to write variant context file to
        variantcontextfile : VariantContextFile
            Variant contexts to write to file
        vcfsamplemap : dict
            Variant files per sample
        bamsamplemap : dict
            Alignment files per sample
        used_donor_vcfs : list of str
            Donor variant files used to create the variant contexts
        used_donor_bams : list of str
            Donor alignment files used to create the variant contexts
        """
        # Write the variant context file itself
        self.vaselogger.info(f"Writing variant contexts to {varcon_outpath}")
        variantcontextfile.write_variant_context_file(varcon_outpath, self.creation_id)

        # Write the variant context statistics file
        self.vaselogger.info(f"Writing variant context statistics to {outpath}varconstats.txt")
        variantcontextfile.write_variant_context_stats(f"{outpath}varconstats.txt", self.creation_id)

        # Write the variant contexts as a BED file
        self.vaselogger.info(f"Writing the variant contexts to a BED file")
        self.write_bed_file(variantcontextfile.get_variant_contexts(), f"{outpath}variantcontexts.bed",
                            self.creation_id)

        # Write a listfile of the used variant (VCF/BCF) files
        self.vaselogger.info(f"Writing the used donor variant files per sample to {outpath}donor_variant_files.txt")
        self.write_used_donor_files(f"{outpath}donorvcfs.txt", vcfsamplemap, used_donor_vcfs, self.creation_id)

        # Write a listfile of the used alignment (BAM/CRAM) files
        self.vaselogger.info(f"Writing the used donor alignment files per sample to {outpath}donor_alignment_files.txt")
        self.write_used_donor_files(f"{outpath}donor_alignment_files.txt", bamsamplemap, used_donor_bams,
                                    self.creation_id)

    def write_pmode_linkfile(self, outpath, context_fqs_links):
        """Writes the link from variant context identifier to fastq output files for P-mode.

        This output file is meant to link a VaSeBuilder run with fastq output files.

        Parameters
        ----------
        outpath : str
            Path to output folder to write the link file to.
        context_fqs_links : dict of str: list
            Fastq files per variant context identifier
        """
        plinkloc = f"{outpath}plink_{self.creation_id}.txt"
        try:
            with open(plinkloc, "w") as plinkfile:
                plinkfile.write(f"#VBUUID: {self.creation_id}\n")
                for contextid, fqfiles in context_fqs_links.items():
                    plinkfile.write(f"{contextid}\t{fqfiles[0]}\t{fqfiles[1]}\n")
        except IOError:
            self.vaselogger.warning(f"Could not write P-mode link file {plinkloc} for run {self.creation_id}")

    def shuffle_donor_read_identifiers(self, donorreads, s=2):
        """Shuffles and returns a list of donor read identifiers.

        Prior to shuffling the donor read identifiers, the provided donor read list is copied to preserve the original.
        The new list is then sorted and a seed is set to ensure that the shuffling can be reproduced if hte same data
        is used.

        Parameters
        ----------
        donorreads :
            Donor read identifiers to shuffle
        s: int
            Seed to set for shuffling (default = 2)

        Returns
        -------
        shuffled_donor_reads : list of str
            Shuffled donor read identifiers
        """
        shuffled_donor_reads = donorreads.copy()
        shuffled_donor_reads.sort()

        random.seed(s)
        random.shuffle(shuffled_donor_reads)
        return shuffled_donor_reads

    def shuffle_donor_add_positions(self, num_of_template_reads, num_of_donor_reads, s=2):
        """Shuffle positions to

        Parameters
        ----------
        num_of_template_reads : int
            Number of template reads in the acceptor
        num_of_donor_reads : int
            Number of donor reads to be added
        s : int
            Seed to set for shuffling (default = 2)

        Returns
        -------
        add_positions : list of int
            Shuffled positions in fastq file to add donor reads to
        """
        random.seed(s)
        self.vaselogger.debug(f"Semi random donor add positions seed set to {s}")
        add_positions = []

        # Check whether the number of donor reads to add exceeds the number of acceptor reads in the template.
        if num_of_donor_reads > num_of_template_reads:
            self.vaselogger.debug("Found more donor reads to add than")
            random_sample_size = num_of_donor_reads
            while random_sample_size > 0:
                pos_to_add = []
                if random_sample_size >= num_of_template_reads:
                    pos_to_add = random.sample(range(0, num_of_template_reads), num_of_template_reads)
                else:
                    pos_to_add = random.sample(range(0, ), random_sample_size)
                add_positions.extend(pos_to_add)
                random_sample_size = random_sample_size - num_of_template_reads
            return add_positions
        add_positions = random.sample(range(0, num_of_template_reads), num_of_donor_reads)
        return add_positions

    def write_donor_output_bam(self, bamoutpath, donorreads, acceptorbam):
        """Writes a set of donor reads as a bam file.

        Parameters
        ----------
        bamoutpath : str
            Path and name to write the output to
        donorreads : list of DonorBamRead
            Donor reads to place in the output BAM file
        acceptorbam
        """
        print("Constructing BAM file")

    def read_is_hard_clipped(self, fetchedread):
        """Returns whether the provided read is hard-clipped

        Parameters
        ----------
        fetchedread : pysam.AlignedSegment
            Read fetched from alignment file with pysam

        Returns
        -------
        bool
            True if read is hardclipped, False if not or read has no cigar string
        """
        if fetchedread.cigarstring is not None:
            return "H" in fetchedread.cigarstring
        return False

    def check_template_size(self, templatefqloc):
        """Returns the number of reads in the template fastq file.

        The returned number is divided by 4 as each read entry consists of four lines.

        Parameters
        ----------
        templatefqloc : str
            Path to acceptor fastq file to check

        Returns
        -------
        int
            Number of read entries in the template fastq file
        """
        line_count = 0
        with gzip.open(templatefqloc, "r") as templatefq:
            for fileline in templatefq:
                line_count += 1
        return int(line_count/4)

    def build_fastq_v2(self, acceptorfq_filepaths, acceptorreads_toskip, donor_context_reads, forward_or_reverse,
                       vasefq_outpath, random_seed, donor_read_insert_data):
        """Builds and writes a set of validation fastq files.

        A set of validation fastq files is build using a set of template/acceptor fastq files. Acceptor reads will be
        filtered from these file via read identifier and donor reads will be added. Donor readsd reads will be added at
        semi random positions. The random.seed() method is used to ensure reproducibility as using the same data (in
        the same order) and the same seed results in the same shuffle.

        Parameters
        ----------
        acceptorfq_filepaths : list of str
            Paths to template fastq files to use
        acceptorreads_toskip : list of str
            Identifiers of acceptor reads to skip
        donor_context_reads : list of tuple
            Donor reads to add to the validation fastq files
        forward_or_reverse : str
            Write forward ('1') or reverse ('2')
        vasefq_outpath :
            Path to write VaSeBuilder validation fastq files to
        """
        # Split all donor reads to add over the template fastq files
        donor_read_ids = [x[0] for x in donor_context_reads]
        donor_read_ids.sort()
        donor_read_ids = list(set(donor_read_ids))
        # self.vaselogger.debug(f"Donor read ids: {donor_read_ids}")
        distributed_read_ids = self.divide_donorfastqs_over_acceptors(donor_read_ids, len(acceptorfq_filepaths))
        # self.vaselogger.debug(f"Distributed read ids: {distributed_read_ids}")

        # Iterate over the R1/R2 fastq in files to use as templates for the
        for x in range(0, len(acceptorfq_filepaths)):
            # Collect the donor reads to write.
            add_donor_ids = distributed_read_ids[x]
            self.vaselogger.debug(f"Distributed read ids: {add_donor_ids}")
            add_donor_reads = [x for x in donor_context_reads if x[0] in add_donor_ids]
            self.vaselogger.debug(f"Will add {len(add_donor_ids)} donor reads")

            # Write the new VaSe FastQ file.
            vasefq_outname = self.set_fastq_out_path(vasefq_outpath, forward_or_reverse, x + 1)
            self.vaselogger.debug(f"Set FastQ output path to: {vasefq_outname}")

            # Check whether to add the fastq file name into the donor read insert position map
            if vasefq_outname.split(".")[0][:-3] not in donor_read_insert_data:
                donor_read_insert_data[vasefq_outname.split(".")[0][:-3]] = {}

            self.write_vase_fastq_v2(acceptorfq_filepaths[x], vasefq_outname,
                                     acceptorreads_toskip, add_donor_reads, add_donor_ids,
                                     forward_or_reverse, random_seed, donor_read_insert_data)

    def write_vase_fastq_v2(self, acceptor_infq, fastq_outpath,
                            acceptorreads_toskip, donorbamreaddata,
                            donor_readids, fr, random_seed, donor_read_insert_data):
        """Creates and writes a single VaSeBuilder validation fastq file.

        Parameters
        ----------
        acceptor_infq : str
            Path to
        fastq_outpath : str
            Path to write VaSeBuilder vaildation fastq file to
        acceptorreads_toskip : list of str
            Acceptor read identifiers to exclude from the validation fastq file
        donorbamreaddata : list of tuple
            Donor reads to writes
        fr : str
            Forward('1') or reverse ('2') fastq file
        """
        fastq_prefix = fastq_outpath.split(".")[0][:-3]
        try:
            fqgz_outfile = io.BufferedWriter(open(fastq_outpath, "wb"))
            self.vaselogger.debug(f"Writing data to validation fastq {fastq_outpath}")

            cur_line_index = 0    # Current read position in the template fastq
            cur_add_index = 0    # Current read position in the validation fastq

            # Determine where to semi randomly add the donor reads in the fastq
            num_of_template_reads = self.check_template_size(acceptor_infq)
            self.vaselogger.debug(f"Template has {num_of_template_reads} reads")
            donor_add_positions = self.shuffle_donor_add_positions(num_of_template_reads, len(donor_readids),
                                                                   random_seed)
            # self.vaselogger.debug(f"Add positions for {fastq_outpath} = {donor_add_positions}")
            donor_reads_to_addpos = self.link_donor_addpos_reads_v2(donor_add_positions, donor_readids,
                                                                    donorbamreaddata)
            # self.vaselogger.debug(f"Read to add pos for {fastq_outpath} = {donor_reads_to_addpos}")

            # Open the template fastq and write filtered data to a new fastq.gz file.
            fqgz_infile = io.BufferedReader(gzip.open(acceptor_infq, "rb"))
            self.vaselogger.debug(f"Opened template FastQ: {acceptor_infq}")
            for fileline in fqgz_infile:

                # Check if we are located at a read identifier.
                if fileline.startswith(b"@"):
                    if fileline.decode("utf-8").split()[0][1:] not in acceptorreads_toskip:
                        fqgz_outfile.write(fileline)
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
                        fqgz_outfile.write(next(fqgz_infile))
                        cur_add_index += 1
                    else:
                        self.vaselogger.debug(f"Skipping acceptor read {fileline}")

                    # Check if we need to add a donor read at the current position
                    if cur_line_index in donor_reads_to_addpos:
                        # self.vaselogger.debug(f"{cur_line_index} is in list of positions to add donor")
                        for donorread in donor_reads_to_addpos[cur_line_index]:
                            if donorread[1] == fr:
                                fqlines = ("@" + str(donorread[0]) + "\n"
                                           + str(donorread[2]) + "\n"
                                           + "+\n"
                                           + str(donorread[3]) + "\n")
                                fqgz_outfile.write(fqlines.encode("utf-8"))
                                cur_add_index += 1
                                self.vaselogger.debug(f"Added donor read {donorread[0]}/{donorread[1]} at "
                                                      f"{cur_add_index}")
                                self.add_donor_insert_data(fastq_prefix, donorread[0], fr, cur_add_index,
                                                           donor_read_insert_data)
                    cur_line_index += 1
            fqgz_infile.close()

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

    def link_donor_addpos_reads_v2(self, donor_addpos, donor_read_ids, donor_reads):
        """Links and returns donor add positions and donor reads.

        Parameters
        ----------
        donor_addpos : list of int
            Positions in validation fastq to add donor reads at
        donor_reads : list of tuple
            Donor reads to add to validation fastq

        Returns
        -------
        add_posread_link : dict of list
            Donor reads to add per add position
        """
        add_posread_link = {}
        for addpos, dread_id in zip(donor_addpos, donor_read_ids):
            if addpos not in add_posread_link:
                add_posread_link[addpos] = []
            add_posread_link[addpos].extend([x for x in donor_reads if x[0] == dread_id])
        return add_posread_link

    def read_donor_fastq(self, donor_fastq, forward_reverse, donor_read_data):
        """Reads and saves the donor fastq reads.

        Parameters
        ----------
        donor_fastq : str
            Path to the donor fastq file to read
        forward_reverse : str
            Reading forward (R1) or reverse (R2) donor read data
        donor_read_data : dict
            Dictionary saving the donor read data

        Returns
        -------
        donor_read_data : dict
            Updated dictionary containing donor read data
        """
        try:
            with gzip.open(donor_fastq, 'rt') as dfqfile:
                for dfqfileline in dfqfile:
                    if dfqfileline.startswith("@"):
                        dfq_readid = dfqfileline[1:].strip()
                        dfq_readseq = next(dfqfile).strip()
                        dfq_optline = next(dfqfile).strip()
                        dfq_readqual = next(dfqfile).strip()
                        if dfq_readid not in donor_read_data:
                            donor_read_data[dfq_readid] = []
                        donor_read_data[dfq_readid].append((dfq_readid, forward_reverse, dfq_readseq, dfq_readqual))
        except IOError:
            self.vaselogger.debug(f"Could not read donor fastq file {donor_fastq}")
        finally:
            return donor_read_data

    def run_ac_mode_v2(self, afq1_in, afq2_in, dfqs, varconfile, random_seed, outpath):
        """Runs VaSeBuilder AC-mode.

        This run mode builds a set of validation fastq files by adding already existing fastq files to acceptor fastq
        files used as template. The variant context file used to construct the variant contexts of the donor fastq files
        is required to filter out the acceptor reads.to be replaced with the donor data.

        Parameters
        ----------
        afq1_in : list of str
            R1 acceptor fastq files to use as template
        afq2_in : list of str
            R2 acceptor fastq files to use as template
        dfqs : list of list of str
            R1 and R2 donor fastq files to add
        varconfile: VariantContextFile
            Variant context to use for filtering out acceptor reads
        outpath: str
            Path to folder to write the output to
        """
        self.vaselogger.info("Running VaSeBuilder AC-mode")
        # Split the donor fastqs into an R1 and R2 group
        r1_dfqs = [dfq[0] for dfq in dfqs]
        r2_dfqs = [dfq[1] for dfq in dfqs]

        # Get the list of acceptor reads to exclude
        skip_list = varconfile.get_all_variant_context_acceptor_read_ids()
        skip_list.sort()

        # Read all the donor read fastq data
        donor_reads = {}
        for fr_i, dfq_i in zip(["1", "2"], [r1_dfqs, r2_dfqs]):
            for x in range(len(dfq_i)):
                donor_reads = self.read_donor_fastq(dfq_i[x], fr_i, donor_reads)
                # self.vaselogger.debug(f"ACmodeV2 Donor reads: {donor_reads}")

        # Divide the donor reads equally over the template fastq files
        donor_read_ids = list(donor_reads.keys())
        donor_read_ids.sort()
        # self.vaselogger.debug(f"ACmodeV2 dread ids: {donor_read_ids}")

        # Distribute the read identifiers over the fastq files
        distributed_donor_read_ids = self.divide_donorfastqs_over_acceptors(donor_read_ids, len(afq1_in))
        # self.vaselogger.debug(f"ACmodeV2 distr dread ids: {distributed_donor_read_ids}")

        # Iterate over the acceptor fastq files
        donor_read_add_data = {}
        for fr_i, afq_i in zip(["1", "2"], [afq1_in, afq2_in]):
            self.build_fastqs_from_donors_v2(afq_i, skip_list, distributed_donor_read_ids, donor_reads, fr_i,
                                             random_seed, outpath, donor_read_add_data)
        self.write_donor_insert_positions_v2(donor_read_add_data, f"{outpath}_donor_read_insert_positions.txt")

    def build_fastqs_from_donors_v2(self, acceptor_fqsin, acceptor_reads_to_exclude, distributed_donor_reads,
                                    donor_reads, forward_reverse, random_seed, outpath, donor_read_insert_data):
        """Builds a set of R1/R2 validation fastq files using existing donor fastq files.

        Parameters
        ----------
        acceptor_fqsin: list of str
            Template/Acceptor fastq files to use
        acceptor_reads_to_exclude: list of str
        distributed_donor_reads : list of list of str
        donor_reads: dict
        forward_reverse: str
            Whether to construct an R1 or R2 fastq file
        random_seed: int
            Random seed to use for shuffling donor add positions
        outpath: str
            Path top write produced fastq file to
        donor_read_insert_data:
        """
        for x in range(len(acceptor_fqsin)):
            fqoutname = self.set_fastq_out_path(outpath, forward_reverse, x + 1)
            self.vaselogger.debug(f"bfmdV2 xdistr dread ids: {distributed_donor_reads[x]}")

            # Add the fastq file entry to the insert positions map
            if fqoutname.split(".")[0][:-3] not in donor_read_insert_data:
                donor_read_insert_data[fqoutname.split(".")[0][:-3]] = {}

            # Filter the donor reads specific to the template file
            selected_donor_reads = [donor_reads[y] for y in distributed_donor_reads[x]]
            donor_reads_to_add = []
            [donor_reads_to_add.extend(x) for x in selected_donor_reads]

            # self.vaselogger.debug(f"bfmdV2 dread ids to add: {donor_reads_to_add}")
            self.write_vase_fastq_v2(acceptor_fqsin[x], fqoutname, acceptor_reads_to_exclude, donor_reads_to_add,
                                     distributed_donor_reads[x], forward_reverse, random_seed, donor_read_insert_data)

    def write_donor_insert_positions_v2(self, inserted_position_data, outpath):
        """Writes the insert positions for each set of reads.

        Insert positions are written per read identifier.

        Parameters
        ----------
        inserted_position_data : dict
            Position data
        outpath : str
            Path and name to write the donor read insert position data to
        """
        try:
            with open(outpath, "w") as ipd_outfile:
                ipd_outfile.write(f"#VBUUID: {self.creation_id}\n")
                ipd_outfile.write("ReadId\tR1_InsertPos\tFastqR1Out\tR2_InsertPos\tFastqR2Out\n")

                # Iterate over the donor read insert position data and write to file
                for fastqout in inserted_position_data:
                    readid_list = [x for x in inserted_position_data[fastqout]]
                    readid_list.sort()

                    # Iterate over each read identifier for the specified output R1 and R2 fastq set.
                    for readid in readid_list:
                        r1_insertpos = self.get_saved_insert_position("1", inserted_position_data[fastqout][readid])
                        r2_insertpos = self.get_saved_insert_position("2", inserted_position_data[fastqout][readid])

                        ipd_outfile.write(f"{readid}\t{r1_insertpos}\t{fastqout}_R1.fastq\t{r2_insertpos}\t"
                                          f"{fastqout}_R2.fastq\n")
        except IOError:
            self.vaselogger.warning(f"Could not write donor insert position data to {outpath}")

    def get_saved_insert_position(self, readpn, read_insert_data):
        """Returns the insert position of a read, or 'NA' if not available.

        Parameters
        ----------
        readpn : str
            The readpair number ('1' or '2')
        read_insert_data : tuple of str and int
            The read insert data

        Returns
        -------
        int or str
            The insert position, 'NA' if not available
        """
        if readpn in read_insert_data:
            readpn_index = read_insert_data.index(readpn)
            if readpn_index < len(read_insert_data)-1:
                return read_insert_data[read_insert_data.index(readpn)+1]
        return "NA"

    def add_donor_insert_data(self, fqoutname, readid, forward_reverse, insertpos, donor_insert_data):
        """

        Parameters
        ----------
        fqoutname : str
            Name of fastq file
        readid : str
            Read identifier for which to save the insert position
        forward_reverse : str
            Indicator whether the read is R1 or R2
        insertpos : int
            Position the donor read was inserted at in the validation fastq
        donor_insert_data: dict
            Saved donor read insertions into validation fastq
        """
        if fqoutname in donor_insert_data:
            if readid in donor_insert_data[fqoutname]:
                if forward_reverse not in donor_insert_data[fqoutname][readid]:
                    donor_insert_data[fqoutname][readid] = donor_insert_data[fqoutname][readid] + \
                                                           (forward_reverse, insertpos)
            else:
                donor_insert_data[fqoutname][readid] = (forward_reverse, insertpos)

    def refetch_donor_reads(self, variant_context_file, donor_alignment_files, genome_reference):
        """Refetches the donor reads from a set of donor BAM files.

        Parameters
        ----------
        variant_context_file : VariantContextFile
            Collection of variant contexts
        donor_alignment_files : dict
            List of paths to donor alignment files
        genome_reference : str
            Path to genome reference file
        """
        varcon_per_sample_id = variant_context_file.get_variant_contexts_by_sampleid()

        # Start iterating over the samples
        for sampleid, donor_aln_file in donor_alignment_files:
            try:
                dalnfile = pysam.AlignmentFile(donor_aln_file, reference_filename=genome_reference)

                # Iterate over the variant contexts for the current sample
                for varcon in varcon_per_sample_id[sampleid]:
                    donor_read_ids = set(varcon.get_donor_read_ids())
                    fetched_reads = self.get_variant_reads(varcon.get_variant_context_id(),
                                                           varcon.get_variant_context_chrom(),
                                                           varcon.get_variant_context_start(),
                                                           varcon.get_variant_context_end(), dalnfile)

                    # Filter out fetched reads not satisfying the donor read identifiers
                    donor_reads = [fdread for fdread in fetched_reads if fdread.query_name in donor_read_ids]
                    variant_context_file.set_variant_context_donor_reads(varcon.get_variant_context_id, donor_reads)
                dalnfile.close()
            except IOError:
                self.vaselogger.warning(f"Could not open donor alignment file {donor_aln_file}")

    def merge_variant_contexts(self, varcon1, varcon2):
        """Merges two variant contexts and returns the mrged context.

        First the associated acceptor and donor contexts are merged. Then the maximum window from both variant contexts
        is determined and reads are merged. The newly combined variant context will retain the context identifier of the
        first variant context.

        Parameters
        ----------
        varcon1 : VariantContext
            First variant context to be merged
        varcon2 : VariantContext
            Second variant context to be merged

        Returns
        -------
        combined_varcon : VariantContext
            New merged variant context
        """
        # Combine the acceptor and donor contexts
        combined_acceptor_context = self.merge_overlap_contexts(varcon1.get_acceptor_context(),
                                                                varcon2.get_acceptor_context())
        combined_donor_context = self.merge_overlap_contexts(varcon1.get_donor_context(),
                                                             varcon2.get_donor_context())

        # Obtain a list of acceptor and donor reads from both variant contexts
        vareads = varcon1.get_acceptor_reads() + varcon2.get_acceptor_reads()
        vdreads = varcon1.get_donor_reads() + varcon2.get_donor_reads()

        # Combine the two variant contexts by determining the new context window and acceptor and donor reads
        combined_window = self.determine_largest_context(varcon1.get_variant_context_origin(), varcon1.get_context(),
                                                         varcon2.get_context())
        combined_vareads = self.uniqify_variant_reads(vareads)
        combined_vdreads = self.uniqify_variant_reads(vdreads)

        # Set the new combined variant context
        combined_varcon = VariantContext(varcon1.get_variant_context_id(), varcon1.get_variant_context_sample(),
                                         *combined_window, combined_vareads, combined_vdreads,
                                         combined_acceptor_context, combined_donor_context)
        return combined_varcon

    def merge_overlap_contexts(self, context1, context2):
        """Merges two acceptor or donor contexts and returns the new merged context.

        Acceptor/Donor contexts are merged by constructing the maximum size of both contexts. Reads in both contexts are
        also merged and uniqified.

        Parameters
        ----------
        context1 : OverlapContext
            First acceptor/donor context to be merged
        context2 : OverlapContext
            Second acceptor/donor context to be merged

        Returns
        -------
        combined_accdon_context : OverlapContext
            Combined acceptor/donor context
        """
        adreads = context1.get_context_bam_reads() + context2.get_context_bam_reads()
        combined_window = self.determine_largest_context(context1.get_context_origin(), context1.get_context(),
                                                         context2.get_context())
        combined_adreads = self.uniqify_variant_reads(adreads)
        combined_accdon_context = OverlapContext(context1.get_context_id(), context1.get_sample_id(), *combined_window,
                                                 combined_adreads)
        return combined_accdon_context
