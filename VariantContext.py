import statistics
from OverlapContext import OverlapContext
from DonorBamRead import DonorBamRead


class VariantContext:
    # Sets the variant context data.
    def __init__(self, varconid, sampleid,
                 varconchrom, varconorigin,
                 varconstart, varconend,
                 acceptorreads, donorreads,
                 acceptor_context=None, donor_context=None):
        self.context_id = varconid
        self.sample_id = sampleid
        self.variant_context_chrom = varconchrom
        self.variant_context_origin = int(varconorigin)
        self.variant_context_start = int(varconstart)
        self.variant_context_end = int(varconend)
        # Saves the acceptor reads that overlap with the entire variant
        # context.
        self.variant_context_areads = acceptorreads
        # Saves the donor reads that overlap with the entire variant
        # context.
        self.variant_context_dreads = donorreads
        # Saves the acceptor context (this is the context from acceptor
        # reads overlapping with the variant itself).
        self.variant_acceptor_context = acceptor_context
        # Saves the donor context (this is the context from donor reads
        # overlapping with the variant itself).
        self.variant_donor_context = donor_context
        self.unmapped_acceptor_mate_ids = []
        self.unmapped_donor_mate_ids = []

    # ===METHODS TO OBTAIN DATA FROM THE VARIANT CONTEXT DATA==================
    # Returns the variant context identifier.
    def get_variant_context_id(self):
        return self.context_id

    # Returns the variant context sample.
    def get_variant_context_sample(self):
        return self.sample_id

    # Returns the variant context chromosome.
    def get_variant_context_chrom(self):
        return self.variant_context_chrom

    # Returns the variant context origin.
    def get_variant_context_origin(self):
        return self.variant_context_origin

    # Returns the variant context start position.
    def get_variant_context_start(self):
        return self.variant_context_start

    # Returns the variant context end position.
    def get_variant_context_end(self):
        return self.variant_context_end

    # Returns the variant context acceptor reads.
    def get_acceptor_reads(self):
        return self.variant_context_areads

    # Returns the variant context donor reads.
    def get_donor_reads(self):
        return self.variant_context_dreads

    # Returns the acceptor context (the context from acceptor reads
    # overlapping the variant itself).
    def get_acceptor_context(self):
        return self.variant_acceptor_context

    # Returns the donor context (the context from donor reads
    # overlapping the variant itself).
    def get_donor_context(self):
        return self.variant_donor_context

    # Returns a list of acceptor reads overlapping with the variant
    # context.
    def get_unmapped_acceptor_mate_ids(self):
        return self.unmapped_acceptor_mate_ids

    # Returns a list of donor reads overlapping with the variant
    # context.
    def get_unmapped_donor_mate_ids(self):
        return self.unmapped_donor_mate_ids

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATING) OF THE===============
    # ===VARIANT CONTEXT===
    # Returns the variant context length.
    def get_variant_context_length(self):
        return abs(self.variant_context_end - self.variant_context_start)

    # Returns the distance of the variant context start position from
    # the variant context origin.
    def get_start_distance_from_origin(self):
        return abs(self.variant_context_origin - self.variant_context_start)

    # Returns the distance of the variant context end position from the
    # variant context origin.
    def get_end_distance_from_origin(self):
        return abs(self.variant_context_end - self.variant_context_origin)

    # ===METHODS TO OBTAIN VARIANT CONTEXT ACCEPTOR READ DATA==================
    # Returns the number of variant context acceptor reads.
    def get_number_of_acceptor_reads(self):
        if self.variant_context_areads is None:
            return 0
        return len(self.variant_context_areads)

    # Returns the identifiers of acceptor reads overlapping with the
    # variant context.
    def get_acceptor_read_ids(self):
        if self.variant_context_areads is None:
            return [None]
        return list(set([x.get_bam_read_id() for x in self.variant_context_areads]))

    # Returns the list of left most acceptor read positions,
    def get_acceptor_read_starts(self):
        if self.variant_context_areads is None:
            return [None]
        return [x.get_bam_read_ref_pos() for x in self.variant_context_areads]

    # Returns a list of the left most positions of the R1 variant
    # context accecptor BAM reads.
    def get_acceptor_read_left_positions(self):
        if self.variant_context_areads is None:
            return [None]
        return [x.get_bam_read_ref_pos()
                for x in self.variant_context_areads if x.is_read1()]

    # Returns the list of all end positions for all variant context
    # acceptor reads.
    def get_acceptor_read_ends(self):
        if self.variant_context_areads is None:
            return [None]
        return [x.get_bam_read_ref_end() for x in self.variant_context_areads]

    # Returns a list of the right most positions ofr the R2 variant
    # context acceptor BAM read.
    def get_acceptor_read_right_positions(self):
        if self.variant_context_areads is None:
            return [None]
        return [x.get_bam_read_ref_end()
                for x in self.variant_context_areads if x.is_read2()]

    # ===METHODS TO OBTAIN VARIANT CONTEXT DONOR READ DATA=====================
    # Returns the number of variant context donor reads.
    def get_number_of_donor_reads(self):
        return len(self.variant_context_dreads)

    # Returns the identifiers of donor reads overlapping with the
    # variant context.
    def get_donor_read_ids(self):
        return list(set([x.get_bam_read_id() for x in self.variant_context_dreads]))

    # Returns the list of variant context donor read starting positions.
    def get_donor_read_starts(self):
        return [x.get_bam_read_ref_pos() for x in self.variant_context_dreads]

    # Returns the list of variant context donor read pairs left most
    # positions (start pos of read 1).
    def get_donor_read_left_positions(self):
        return [x.get_bam_read_ref_pos()
                for x in self.variant_context_dreads if (x.is_read1())]

    # Returns a list of all donor read ending positions
    def get_donor_read_ends(self):
        return [x.get_bam_read_ref_end() for x in self.variant_context_dreads]

    # Returns a list of all variant context donor reads right most
    # positions (end pos of read 2).
    def get_donor_read_right_positions(self):
        return [x.get_bam_read_ref_end()
                for x in self.variant_context_dreads if (x.is_read2())]

    # ===METHODS TO ADD DATA TO THE VARIANT CONTEXT============================
    # Sets the acceptor context of the variant context.
    def set_acceptor_context(self, acceptor_context):
        self.variant_acceptor_context = acceptor_context

    # Sets the donor context of the variant context.
    def set_donor_context(self, donor_context):
        self.variant_donor_context = donor_context

    # Adds an acceptor context to the variant context by creating it
    # from data values.
    def add_acceptor_context(self, contextid, sampleid,
                             contextchrom, contextorigin,
                             contextstart, contextend,
                             acceptorreads):
        self.variant_acceptor_context = OverlapContext(
                contextid, sampleid,
                contextchrom, contextorigin,
                contextstart, contextend,
                acceptorreads
                )

    # Adds a donor context to the variant context by creating it from
    # data values.
    def add_donor_context(self, contextid, sampleid,
                          contextchrom, contextorigin,
                          contextstart, contextend,
                          donorreads):
        self.variant_donor_context = OverlapContext(
                contextid, sampleid,
                contextchrom, contextorigin,
                contextstart, contextend,
                donorreads
                )

    # ===METHODS TO OBTAIN VARIANT CONTEXT UNMAPPED MATE READ DATA=============
    # Returns the variant context acceptor read ids that have an
    # unmapped mate.
    def get_unmapped_acceptor_read_ids(self):
        return self.unmapped_acceptor_mate_ids

    # Returns the variant context donor read ids that have an unmapped
    # mate.
    def get_unmapped_donor_read_ids(self):
        return self.unmapped_donor_mate_ids

    # Adds a variant context appector mate identifier.
    def add_unmapped_acceptor_mate_id(self, mateid):
        self.unmapped_acceptor_mate_ids.append(mateid)

    # Adds a variant context donor mate identifier.
    def add_unmapped_donor_mate_id(self, mateid):
        self.unmapped_donor_mate_ids.append(mateid)

    # Sets the variant context unmapped acceptor mate ids.
    def set_unmapped_acceptor_mate_ids(self, mateids):
        self.unmapped_acceptor_mate_ids = mateids

    # Sets the variant context unmapped donor mate ids.
    def set_unmapped_donor_mate_ids(self, mateids):
        self.unmapped_donor_mate_ids = mateids

    # Returns whether a specified variant context acceptor read has an
    # unmapped mate.
    def acceptor_read_has_unmapped_mate(self, readid):
        return readid in self.unmapped_acceptor_mate_ids

    # Returns whether a specified variant context donor read has an
    # unmapped mate.
    def donor_read_has_unmapped_mate(self, readid):
        return readid in self.unmapped_donor_mate_ids

    # Returns the number of variant context acceptor reads with an
    # unmapped mate.
    def get_number_of_unmapped_acceptor_mates(self):
        return len(self.unmapped_acceptor_mate_ids)

    # Returns the number of variant context donor reads with an
    # unmapped mate.
    def get_number_of_unmapped_donor_mates(self):
        return len(self.unmapped_donor_mate_ids)

    # ===METHODS TO ADD UNMAPPED MATES TO THE ACCEPTOR AND DONOR CONTEXT=======
    # Sets the unmapped mate ids for the acceptor context.
    def set_acceptor_context_unmapped_mates(self, mateids):
        self.variant_acceptor_context.set_unmapped_mate_ids(mateids)

    # Adds an unmapped read id to the acceptor context.
    def add_acceptor_context_unmapped_mate(self, ureadid):
        self.variant_acceptor_context.add_unmapped_mate_id(ureadid)

    # Sets the unmapped mate ids for the donor context.
    def set_donor_context_unmapped_mates(self, mateids):
        self.variant_donor_context.set_unmapped_mate_ids(mateids)

    # Adds an unmapped read id to the donor context.
    def add_donor_context_unmapped_mate_id(self, ureadid):
        self.variant_donor_context.setUnmappedMateId(ureadid)

    # ===METHODS TO OBTAIN STATISTICS OF THE VARIANT CONTEXT===================
    # Returns the average and median quality of the acceptor reads
    # associated with
    def get_average_and_median_acceptor_read_qual(self):
        return self.get_average_and_median_read_qual(self.variant_context_areads)

    # Returns the average and median quality of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_qual(self):
        return self.get_average_and_median_read_qual(self.variant_context_dreads)

    # Returns the average and median read quality.
    def get_average_and_median_read_qual(self, contextreads):
        if contextreads is not None:
            avgmedqual = []
            for contextread in contextreads:
                avgmedqual.append(contextread.get_average_qscore())
            return ([statistics.mean(avgmedqual),
                     statistics.median(avgmedqual)])
        return [None, None]

    # Returns the average and median mapq values of the acceptor reads
    # associated with this variant context.
    def get_average_and_median_acceptor_read_mapq(self):
        return self.get_average_and_median_read_mapq(self.variant_context_areads)

    # Returns the average and median mapq values of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_mapq(self):
        return self.get_average_and_median_read_mapq(self.variant_context_dreads)

    # Returns the average and median read MapQ of this variant context.
    def get_average_and_median_read_mapq(self, contextreads):
        if contextreads is not None:
            avgmedmapq = []
            for contextread in contextreads:
                avgmedmapq.append(contextread.get_mapping_qual())
            return ([statistics.mean(avgmedmapq),
                     statistics.median(avgmedmapq)])
        return [None, None]

    # Returns the average and median length of the acceptor reads
    # associated with this variant context.
    def get_average_and_median_acceptor_read_length(self):
        return self.get_average_and_median_read_length(self.variant_context_areads)

    # Returns the average and median length of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_length(self):
        return self.get_average_and_median_read_length(self.variant_context_dreads)

    # Returns the average and median read length.
    def get_average_and_median_read_length(self, contextreads):
        if contextreads is not None:
            avgmedlen = []
            for contextread in contextreads:
                if contextread.get_bam_read_length() is not None:
                    avgmedlen.append(contextread.get_bam_read_length())
            return [statistics.mean(avgmedlen), statistics.median(avgmedlen)]
        return [None, None]

    # ===METHODS TO OBTAIN ACCEPTOR CONTEXT DATA===============================
    # Returns whether the variant context has an acceptor context
    def has_acceptor_context(self):
        return self.variant_acceptor_context is not None

    # Returns the acceptor context identifier (should be the same as the
    # variant context id).
    def get_acceptor_context_id(self):
        return self.variant_acceptor_context.get_context_id()

    # Returns the acceptor context sample id (should be the same as the
    # variant context sample id).
    def get_acceptor_context_sample_id(self):
        return self.variant_acceptor_context.get_sample_id()

    # Returns the chromosome of the acceptor context.
    def get_acceptor_context_chrom(self):
        return self.variant_acceptor_context.get_context_chrom()

    # Returns the origin position of the acceptor context.
    def get_acceptor_context_origin(self):
        return self.variant_acceptor_context.get_context_origin()

    # Returns the starting position of the acceptor context.
    def get_acceptor_context_start(self):
        return self.variant_acceptor_context.get_context_start()

    # Returns the ending position of the acceptor context.
    def get_acceptor_context_end(self):
        return self.variant_acceptor_context.get_context_end()

    # Returns the length of the acceptor context.
    def get_acceptor_context_length(self):
        return self.variant_acceptor_context.get_context_length()

    # Returns the acceptor context reads.
    def get_acceptor_context_reads(self):
        return self.variant_acceptor_context.get_context_bam_reads()

    # Returns the lost of read ids in the acceptor context.
    def get_acceptor_context_read_ids(self):
        return self.variant_acceptor_context.get_context_bam_read_ids()

    # Returns a list of all acceptor context BAM read start positions.
    def get_acceptor_context_read_starts(self):
        return self.variant_acceptor_context.get_context_bam_read_starts()

    # Returns a list of all acceptor context R1 BAM read left positons.
    def get_acceptor_context_read_left_positions(self):
        return self.variant_acceptor_context.get_context_bam_read_left_positions()

    # Returns a list of all acceptor context BAM read end positions.
    def get_acceptor_context_read_ends(self):
        return self.variant_acceptor_context.get_context_bam_read_ends()

    # Returns a list of all acceptor context R2 BAM read end positions.
    def get_acceptor_context_read_right_positions(self):
        return self.variant_acceptor_context.get_context_bam_read_right_positions()

    # Returns a list of all acceptor context BAM read lengths.
    def get_acceptor_context_read_lengths(self):
        return self.variant_acceptor_context.get_context_bam_read_lengths()

    # Returns the list of acceptor context unmapped mate read ids
    def get_acceptor_context_unmapped_mate_ids(self):
        return self.variant_acceptor_context.get_unmapped_read_mate_ids()

    # ===METHODS TO OBTAIN DONOR CONTEXT DATA==================================
    # Returns whether the variant context has a donor context
    def has_donor_context(self):
        return self.variant_donor_context is not None

    # Returns the donor context identifier (should be the same as the
    # variant context id).
    def get_donor_context_id(self):
        return self.variant_donor_context.get_context_id()

    # Returns the acceptor context sample id (should be the same as the
    # variant context sample id).
    def get_donor_context_sample_id(self):
        return self.variant_donor_context.get_sample_id()

    # Returns the chromosome of the donor context.
    def get_donor_context_chrom(self):
        return self.variant_donor_context.get_context_chrom()

    # Returns the origin position of the donor context.
    def get_donor_context_origin(self):
        return self.variant_donor_context.get_context_origin()

    # Returns the starting position of the donor context.
    def get_donor_context_start(self):
        return self.variant_donor_context.get_context_start()

    # Returns the ending position of the donor context.
    def get_donor_context_end(self):
        return self.variant_donor_context.get_context_end()

    # Returns the length of the acceptor context.
    def get_donor_context_length(self):
        return self.variant_donor_context.get_context_length()

    # Returns the donor context reads.
    def get_donor_context_reads(self):
        return self.variant_donor_context.get_context_bam_reads()

    # Returns the donor context read identifers.
    def get_donor_context_read_ids(self):
        return self.variant_donor_context.get_context_bam_read_ids()

    # Returns a list of donor context read start positions.
    def get_donor_context_read_starts(self):
        return self.variant_donor_context.get_context_bam_read_starts()

    # Returns a list of all acceptor context R1 BAM read left positons.
    def get_donor_context_read_left_positions(self):
        return self.variant_donor_context.get_context_bam_read_left_positions()

    # Returns a list of donor context read end positions.
    def get_donor_context_read_ends(self):
        return self.variant_donor_context.get_context_bam_read_ends()

    # Returns a list of all acceptor context R2 BAM read right positons.
    def get_donor_context_read_right_positions(self):
        return self.variant_donor_context.get_context_bam_read_right_positions()

    # Returns a list of donor context BAM read lengths.
    def get_donor_context_read_lengths(self):
        return self.variant_donor_context.get_context_bam_read_lengths()

    # Returns the list of donor context unmapped mate read ids.
    def get_donor_context_unmapped_mate_ids(self):
        return self.variant_donor_context.get_unmapped_read_mate_ids()

    # ===METHODS TO PRODUCE SOME OUTPUT ABOUT THE VARIANT CONTEXT==============
    # Returns a varcon.txt string representation of the variant context.
    def to_string(self):
        if self.variant_context_areads is None:
            ad_ratio = "N/A"
            list_areads = None
            acon_len = None
            aread_count = 0
        else:
            ad_ratio = float(len(self.variant_context_areads)
                             / len(self.variant_context_dreads))
            list_areads = ";".join(list(set(self.get_acceptor_read_ids())))
            acon_len = self.variant_acceptor_context.get_context_length()
            aread_count = len(self.variant_context_areads)
        return (str(self.context_id) + "\t"
                + str(self.sample_id) + "\t"
                + str(self.variant_context_chrom) + "\t"
                + str(self.variant_context_origin) + "\t"
                + str(self.variant_context_start) + "\t"
                + str(self.variant_context_end) + "\t"
                + str(acon_len) + "\t"
                + str(self.variant_donor_context.get_context_length()) + "\t"
                + str(aread_count) + "\t"
                + str(len(self.variant_context_dreads)) + "\t"
                + str(ad_ratio) + "\t"
                + str(list_areads) + "\t"
                + ";".join(list(set(self.get_donor_read_ids()))))

    # Returns a varconstats.txt string representation of the variant
    # context.
    def to_statistics_string(self):
        areadlen = self.get_average_and_median_acceptor_read_length()
        dreadlen = self.get_average_and_median_donor_read_length()
        areadqual = self.get_average_and_median_acceptor_read_qual()
        dreadqual = self.get_average_and_median_donor_read_qual()
        areadmapq = self.get_average_and_median_acceptor_read_mapq()
        dreadmapq = self.get_average_and_median_donor_read_mapq()
        return (str(self.context_id) + "\t"
                + str(areadlen[0]) + "\t"
                + str(dreadlen[0]) + "\t"
                + str(areadlen[1]) + "\t"
                + str(dreadlen[1]) + "\t"
                + str(areadqual[0]) + "\t"
                + str(dreadqual[0]) + "\t"
                + str(areadqual[1]) + "\t"
                + str(dreadqual[1]) + "\t"
                + str(areadmapq[0]) + "\t"
                + str(dreadmapq[0]) + "\t"
                + str(areadmapq[1]) + "\t"
                + str(dreadmapq[1]))
