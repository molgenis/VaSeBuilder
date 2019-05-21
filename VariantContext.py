import statistics
from OverlapContext import OverlapContext
from DonorBamRead import DonorBamRead


class VariantContext:
    # Sets the variant context data.
    def __init__(self, varconId, varconSample,
                 varconChrom, varconOrigin,
                 varconStart, varconEnd,
                 acceptorReads, donorReads,
                 acceptorContext=None, donorContext=None):
        self.contextId = varconId
        self.sampleId = varconSample
        self.variantContextChrom = varconChrom
        self.variantContextOrigin = varconOrigin
        self.variantContextStart = varconStart
        self.variantContextEnd = varconEnd
        # Saves the acceptor reads that overlap with the entire variant
        # context.
        self.variantContextAReads = acceptorReads
        # Saves the donor reads that overlap with the entire variant
        # context.
        self.variantContextDReads = donorReads
        # Saves the acceptor context (this is the context from acceptor
        # reads overlapping with the variant itself).
        self.variantAcceptorContext = acceptorContext
        # Saves the donor context (this is the context from donor reads
        # overlapping with the variant itself).
        self.variantDonorContext = donorContext
        self.unmappedAcceptorMateIds = []
        self.unmappedDonorMateIds = []

    # ===METHODS TO OBTAIN DATA FROM THE VARIANT CONTEXT DATA==================
    # Returns the variant context identifier.
    def get_variant_context_id(self):
        return self.contextId

    # Returns the variant context sample.
    def get_variant_context_sample(self):
        return self.sampleId

    # Returns the variant context chromosome.
    def get_variant_context_chrom(self):
        return self.variantContextChrom

    # Returns the variant context origin.
    def get_variant_context_origin(self):
        return self.variantContextOrigin

    # Returns the variant context start position.
    def get_variant_context_start(self):
        return self.variantContextStart

    # Returns the variant context end position.
    def get_variant_context_end(self):
        return self.variantContextEnd

    # Returns the variant context acceptor reads.
    def get_acceptor_reads(self):
        return self.variantContextAReads

    # Returns the variant context donor reads.
    def get_donor_reads(self):
        return self.variantContextDReads

    # Returns the acceptor context (the context from acceptor reads
    # overlapping the variant itself).
    def get_acceptor_context(self):
        return self.variantAcceptorContext

    # Returns the donor context (the context from donor reads
    # overlapping the variant itself).
    def get_donor_context(self):
        return self.variantDonorContext

    # Returns a list of acceptor reads overlapping with the variant
    # context.
    def get_unmapped_acceptor_mate_ids(self):
        return self.unmappedAcceptorMateIds

    # Returns a list of donor reads overlapping with the variant
    # context.
    def get_unmapped_donor_mate_ids(self):
        return self.unmappedDonorMateIds

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATING) OF THE===============
    # ===VARIANT CONTEXT===
    # Returns the variant context length.
    def get_variant_context_length(self):
        return abs(self.variantContextEnd - self.variantContextStart)

    # Returns the distance of the variant context start position from
    # the variant context origin.
    def get_start_distance_from_origin(self):
        return abs(self.variantContextOrigin - self.variantContextStart)

    # Returns the distance of the variant context end position from the
    # variant context origin.
    def get_end_distance_from_origin(self):
        return abs(self.variantContextEnd - self.variantContextOrigin)

    # ===METHODS TO OBTAIN VARIANT CONTEXT ACCEPTOR READ DATA==================
    # Returns the number of variant context acceptor reads.
    def get_number_of_acceptor_reads(self):
        return len(self.variantContextAReads)

    # Returns the identifiers of acceptor reads overlapping with the
    # variant context.
    def get_acceptor_read_ids(self):
        return [x.get_bam_read_id() for x in self.variantContextAReads]

    # Returns the list of left most acceptor read positions,
    def get_acceptor_read_starts(self):
        return [x.get_bam_read_ref_pos() for x in self.variantContextAReads]

    # Returns a list of the left most positions of the R1 variant
    # context accecptor BAM reads.
    def get_acceptor_read_left_positions(self):
        return [x.get_bam_read_ref_pos()
                for x in self.variantContextAReads if x.is_read1()]

    # Returns the list of all end positions for all variant context
    # acceptor reads.
    def get_acceptor_read_ends(self):
        return [x.get_bam_read_ref_end() for x in self.variantContextAReads]

    # Returns a list of the right most positions ofr the R2 variant
    # context acceptor BAM read.
    def get_acceptor_read_right_positions(self):
        return [x.get_bam_read_ref_end()
                for x in self.variantContextAReads if x.is_read2()]

    # ===METHODS TO OBTAIN VARIANT CONTEXT DONOR READ DATA=====================
    # Returns the number of variant context donor reads.
    def get_number_of_donor_reads(self):
        return len(self.variantContextDReads)

    # Returns the identifiers of donor reads overlapping with the
    # variant context.
    def get_donor_read_ids(self):
        return [x.get_bam_read_id() for x in self.variantContextDReads]

    # Returns the list of variant context donor read starting positions.
    def get_donor_read_starts(self):
        return [x.get_bam_read_ref_pos() for x in self.variantContextDReads]

    # Returns the list of variant context donor read pairs left most
    # positions (start pos of read 1).
    def get_donor_read_left_positions(self):
        return [x.get_bam_read_ref_pos()
                for x in self.variantContextDReads if (x.is_read1())]

    # Returns a list of all donor read ending positions
    def get_donor_read_ends(self):
        return [x.get_bam_read_ref_end() for x in self.variantContextDReads]

    # Returns a list of all variant context donor reads right most
    # positions (end pos of read 2).
    def get_donor_read_right_positions(self):
        return [x.get_bam_read_ref_end()
                for x in self.variantContextDReads if (x.is_read2())]

    # ===METHODS TO ADD DATA TO THE VARIANT CONTEXT============================
    # Sets the acceptor context of the variant context.
    def set_acceptor_context(self, acceptorContext):
        self.variantAcceptorContext = acceptorContext

    # Sets the donor context of the variant context.
    def set_donor_context(self, donorContext):
        self.variantDonorContext = donorContext

    # Adds an acceptor context to the variant context by creating it
    # from data values.
    def add_acceptor_context(self, contextId, sampleId,
                             contextChrom, contextOrigin,
                             contextStart, contextEnd,
                             acceptorReads):
        self.variantAcceptorContext = OverlapContext(
                contextId, sampleId,
                contextChrom, contextOrigin,
                contextStart, contextEnd,
                acceptorReads
                )

    # Adds a donor context to the variant context by creating it from
    # data values.
    def add_donor_context(self, contextId, sampleId,
                          contextChrom, contextOrigin,
                          contextStart, contextEnd,
                          donorReads):
        self.variantDonorContext = OverlapContext(
                contextId, sampleId,
                contextChrom, contextOrigin,
                contextStart, contextEnd,
                donorReads
                )

    # ===METHODS TO OBTAIN VARIANT CONTEXT UNMAPPED MATE READ DATA=============
    # Returns the variant context acceptor read ids that have an
    # unmapped mate.
    def get_unmapped_acceptor_read_ids(self):
        return self.unmappedAcceptorMateIds

    # Returns the variant context donor read ids that have an unmapped
    # mate.
    def get_unmapped_donor_read_ids(self):
        return self.unmappedDonorMateIds

    # Adds a variant context appector mate identifier.
    def add_unmapped_acceptor_mate_id(self, mateId):
        self.unmappedAcceptorMateIds.append(mateId)

    # Adds a variant context donor mate identifier.
    def add_unmapped_donor_mate_id(self, mateId):
        self.unmappedDonorMateIds.append(mateId)

    # Sets the variant context unmapped acceptor mate ids.
    def set_unmapped_acceptor_mate_ids(self, mateIds):
        self.unmappedAcceptorMateIds = mateIds

    # Sets the variant context unmapped donor mate ids.
    def set_unmapped_donor_mate_ids(self, mateIds):
        self.unmappedDonorMateIds = mateIds

    # Returns whether a specified variant context acceptor read has an
    # unmapped mate.
    def acceptor_read_has_unmapped_mate(self, readId):
        return readId in self.unmappedAcceptorMateIds

    # Returns whether a specified variant context donor read has an
    # unmapped mate.
    def donor_read_has_unmapped_mate(self, readId):
        return readId in self.unmappedDonorMateIds

    # Returns the number of variant context acceptor reads with an
    # unmapped mate.
    def get_number_of_unmapped_acceptor_mates(self):
        return len(self.unmappedAcceptorMateIds)

    # Returns the number of variant context donor reads with an
    # unmapped mate.
    def get_number_of_unmapped_donor_mates(self):
        return len(self.unmappedDonorMateIds)

    # ===METHODS TO ADD UNMAPPED MATES TO THE ACCEPTOR AND DONOR CONTEXT=======
    # Sets the unmapped mate ids for the acceptor context.
    def set_acceptor_context_unmapped_mates(self, mateIds):
        self.variantAcceptorContext.set_unmapped_mate_ids(mateIds)

    # Adds an unmapped read id to the acceptor context.
    def add_acceptor_context_unmapped_mate(self, uReadId):
        self.variantAcceptorContext.add_unmapped_mate_id(uReadId)

    # Sets the unmapped mate ids for the donor context.
    def set_donor_context_unmapped_mates(self, mateIds):
        self.variantDonorContext.set_unmapped_mate_ids(mateIds)

    # Adds an unmapped read id to the donor context.
    def add_donor_context_unmapped_mate_id(self, uReadId):
        self.variantDonorContext.setUnmappedMateId(uReadId)

    # ===METHODS TO OBTAIN STATISTICS OF THE VARIANT CONTEXT===================
    # Returns the average and median quality of the acceptor reads
    # associated with
    def get_average_and_median_acceptor_read_qual(self):
        return self.get_average_and_median_read_qual(self.variantContextAReads)

    # Returns the average and median quality of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_qual(self):
        return self.get_average_and_median_read_qual(self.variantContextDReads)

    # Returns the average and median read quality.
    def get_average_and_median_read_qual(self, contextReads):
        if (contextReads is not None):
            avgMedQual = []
            for contextread in contextReads:
                avgMedQual.append(contextread.get_average_qscore())
            return ([statistics.mean(avgMedQual),
                     statistics.median(avgMedQual)])
        return None

    # Returns the average and median mapq values of the acceptor reads
    # associated with this variant context.
    def get_average_and_median_acceptor_read_mapq(self):
        return self.get_average_and_median_read_mapq(self.variantContextAReads)

    # Returns the average and median mapq values of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_mapq(self):
        return self.get_average_and_median_read_mapq(self.variantContextDReads)

    # Returns the average and median read MapQ of this variant context.
    def get_average_and_median_read_mapq(self, contextReads):
        if (contextReads is not None):
            avgMedMapQ = []
            for contextread in contextReads:
                avgMedMapQ.append(contextread.get_mapping_qual())
            return ([statistics.mean(avgMedMapQ),
                     statistics.median(avgMedMapQ)])
        return None

    # Returns the average and median length of the acceptor reads
    # associated with this variant context.
    def get_average_and_median_acceptor_read_length(self):
        return self.get_average_and_median_read_length(self.variantContextAReads)

    # Returns the average and median length of the donor reads
    # associated with this variant context.
    def get_average_and_median_donor_read_length(self):
        return self.get_average_and_median_read_length(self.variantContextDReads)

    # Returns the average and median read length.
    def get_average_and_median_read_length(self, contextReads):
        if (contextReads is not None):
            avgMedLen = []
            for contextread in contextReads:
                avgMedLen.append(contextread.get_bam_read_length())
            return ([statistics.mean(avgMedLen), statistics.median(avgMedLen)])
        return None

    # ===METHODS TO OBTAIN ACCEPTOR CONTEXT DATA===============================
    # Returns the acceptor context identifier (should be the same as the
    # variant context id).
    def get_acceptor_context_id(self):
        return self.variantAcceptorContext.get_context_id()

    # Returns the acceptor context sample id (should be the same as the
    # variant context sample id).
    def get_acceptor_context_sample_id(self):
        return self.variantAcceptorContext.get_sample_id()

    # Returns the chromosome of the acceptor context.
    def get_acceptor_context_chrom(self):
        return self.variantAcceptorContext.get_context_chrom()

    # Returns the origin position of the acceptor context.
    def get_acceptor_context_origin(self):
        return self.variantAcceptorContext.get_context_origin()

    # Returns the starting position of the acceptor context.
    def get_acceptor_context_start(self):
        return self.variantAcceptorContext.get_context_start()

    # Returns the ending position of the acceptor context.
    def get_acceptor_context_end(self):
        return self.variantAcceptorContext.get_context_end()

    # Returns the length of the acceptor context.
    def get_acceptor_context_length(self):
        return self.variantAcceptorContext.get_context_length()

    # Returns the acceptor context reads.
    def get_acceptor_context_reads(self):
        return self.variantAcceptorContext.get_context_bam_reads()

    # Returns the lost of read ids in the acceptor context.
    def get_acceptor_context_read_ids(self):
        return self.variantAcceptorContext.get_context_bam_read_ids()

    # Returns a list of all acceptor context BAM read start positions.
    def get_acceptor_context_read_starts(self):
        return self.variantAcceptorContext.get_context_bam_read_starts()

    # Returns a list of all acceptor context R1 BAM read left positons.
    def get_acceptor_context_read_left_positions(self):
        return self.variantAcceptorContext.get_context_bam_read_left_positions()

    # Returns a list of all acceptor context BAM read end positions.
    def get_acceptor_context_read_ends(self):
        return self.variantAcceptorContext.get_context_bam_read_ends()

    # Returns a list of all acceptor context R2 BAM read end positions.
    def get_acceptor_context_read_right_positions(self):
        return self.variantAcceptorContext.get_context_bam_read_right_positions()

    # Returns a list of all acceptor context BAM read lengths.
    def get_acceptor_context_read_lengths(self):
        return self.variantAcceptorContext.get_context_bam_read_lengths()
    
    # Returns the list of acceptor context unmapped mate read ids
    def get_acceptor_context_unmapped_mate_ids(self):
        return self.variantAcceptorContext.get_unmapped_read_mate_ids()

    # ===METHODS TO OBTAIN DONOR CONTEXT DATA==================================
    # Returns the donor context identifier (should be the same as the
    # variant context id).
    def get_donor_context_id(self):
        return self.variantDonorContext.get_context_id()

    # Returns the acceptor context sample id (should be the same as the
    # variant context sample id).
    def get_donor_context_sample_id(self):
        return self.variantDonorContext.get_sample_id()

    # Returns the chromosome of the donor context.
    def get_donor_context_chrom(self):
        return self.variantDonorContext.get_context_chrom()

    # Returns the origin position of the donor context.
    def get_donor_context_origin(self):
        return self.variantDonorContext.get_context_origin()

    # Returns the starting position of the donor context.
    def get_donor_context_start(self):
        return self.variantDonorContext.get_context_start()

    # Returns the ending position of the donor context.
    def get_donor_context_end(self):
        return self.variantDonorContext.get_context_end()

    # Returns the length of the acceptor context.
    def get_donor_context_length(self):
        return self.variantDonorContext.get_context_length()

    # Returns the donor context reads.
    def get_donor_context_reads(self):
        return self.variantDonorContext.get_context_bam_reads()

    # Returns the donor context read identifers.
    def get_donor_context_read_ids(self):
        return self.variantDonorContext.get_context_bam_read_ids()

    # Returns a list of donor context read start positions.
    def get_donor_context_read_starts(self):
        return self.variantDonorContext.get_context_bam_read_starts()

    # Returns a list of all acceptor context R1 BAM read left positons.
    def get_donor_context_read_left_positions(self):
        return self.variantDonorContext.get_context_bam_read_left_positions()

    # Returns a list of donor context read end positions.
    def get_donor_context_read_ends(self):
        return self.variantDonorContext.get_context_bam_read_ends()

    # Returns a list of all acceptor context R2 BAM read right positons.
    def get_donor_context_read_right_positions(self):
        return self.variantDonorContext.get_context_bam_read_right_positions()

    # Returns a list of donor context BAM read lengths.
    def get_donor_context_read_lengths(self):
        return self.variantDonorContext.get_context_bam_read_lengths()

    # Returns the list of donor context unmapped mate read ids.
    def get_donor_context_unmapped_mate_ids(self):
        return self.variantDonorContext.get_unmapped_read_mate_ids()

    # ===METHODS TO PRODUCE SOME OUTPUT ABOUT THE VARIANT CONTEXT==============
    # Returns a varcon.txt string representation of the variant context.
    def to_string(self):
        return (str(self.contextId) + "\t"
                + str(self.sampleId) + "\t"
                + str(self.variantContextChrom) + "\t"
                + str(self.variantContextOrigin) + "\t"
                + str(self.variantContextStart) + "\t"
                + str(self.variantContextEnd) + "\t"
                + str(self.variantAcceptorContext.get_context_length()) + "\t"
                + str(self.variantDonorContext.get_context_length()) + "\t"
                + str(len(self.variantContextAReads)) + "\t"
                + str(len(self.variantContextDReads)) + "\t"
                + str(float(len(self.variantContextAReads)
                            / len(self.variantContextDReads))) + "\t"
                + ';'.join(self.get_acceptor_read_ids()) + "\t"
                + ';'.join(self.get_donor_read_ids()))

    # Returns a varconstats.txt string representation of the variant
    # context.
    def to_statistics_string(self):
        aReadLen = self.get_average_and_median_acceptor_read_length()
        dReadLen = self.get_average_and_median_donor_read_length()
        aReadQual = self.get_average_and_median_acceptor_read_qual()
        dReadQual = self.get_average_and_median_donor_read_qual()
        aReadMapQ = self.get_average_and_median_acceptor_read_mapq()
        dReadMapQ = self.get_average_and_median_donor_read_mapq()
        return (str(self.contextId) + "\t"
                + str(aReadLen[0]) + "\t"
                + str(dReadLen[0]) + ""
                + str(aReadLen[1]) + ""
                + str(dReadLen[1]) + "\t"
                + str(aReadQual[0]) + "\t"
                + str(dReadQual[0]) + "\t"
                + str(aReadQual[1]) + "\t"
                + str(dReadQual[1]) + "\t"
                + str(aReadMapQ[0]) + "\t"
                + str(dReadMapQ[0]) + "\t"
                + str(aReadMapQ[1]) + "\t"
                + str(dReadMapQ[1]))
