import statistics
from DonorBamRead import DonorBamRead


class OverlapContext:
    # Saves the data associated with the context overlapping with the
    # variant.
    def __init__(self, variantid, sampleid, ovconchrom, ovconorigin,
                 ovconstart, ovconend, bamreads):
        self.contextId = variantid
        self.sampleId = sampleid
        self.contextChrom = ovconchrom
        self.contextOrigin = ovconorigin
        self.contextStart = ovconstart
        self.contextEnd = ovconend
        self.contextBamReads = bamreads
        self.unmappedReadMateIds = []

    # ===METHODS TO GET DATA OF THE OVERLAP CONTEXT============================
    # Returns the context identifier.
    def get_context_id(self):
        return self.contextId

    # Returns the sample identifier the context was based on.
    def get_sample_id(self):
        return self.sampleId

    # Returns the context chromosome name the context is located on.
    def get_context_chrom(self):
        return self.contextChrom

    # Returns the origin position of the context.
    def get_context_origin(self):
        return self.contextOrigin

    # Returns the context start position.
    def get_context_start(self):
        return self.contextStart

    # Returns the end position of the context.
    def get_context_end(self):
        return self.contextEnd

    # Returns the bam reads associated with the context as a list of
    # BamRead objects.
    def get_context_bam_reads(self):
        return self.contextBamReads

    # Returns the list of BAM read ids that have an unmapped mate.
    def get_unmapped_read_mate_ids(self):
        return self.unmappedReadMateIds

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATION) FROM THE=============
    # ===OVERLAP CONTEXT===
    # Returns the lengt of the context.
    def get_context_length(self):
        return abs(self.contextEnd - self.contextStart)

    # Returns the distance of the context start from the context origin.
    def get_start_distance_from_origin(self):
        return abs(self.contextOrigin - self.contextStart)

    # Returns the distance of the context end from the context origin.
    def get_end_distance_from_origin(self):
        return abs(self.contextEnd - self.contextOrigin)

    # ===METHODS TO OBTAIN CONTEXT READ INFORMATION============================
    # Returns the number of saved context reads.
    def get_number_of_context_reads(self):
        return len(self.contextBamReads)

    # Returns a list of BAM read identifiers in the current context.
    def get_context_bam_read_ids(self):
        return [x.get_bam_read_id() for x in self.contextBamReads]

    # Returns a list of all left positions for all BAM reads.
    def get_context_bam_read_starts(self):
        return [x.get_bam_read_ref_pos() for x in self.contextBamReads]

    # Returns a list of all left positions for all R1 BAM reads.
    def get_context_bam_read_left_positions(self):
        return [x.get_bam_read_ref_pos()
                for x in self.contextBamReads if (x.is_read1())]

    # Returns a list of BAM read ending positions for all BAM reads.
    def get_context_bam_read_ends(self):
        return [x.get_bam_read_ref_end() for x in self.contextBamReads]

    # Returns a list of all right positions for all R2 BAM reads.
    def get_context_bam_read_right_positions(self):
        return [x.get_bam_read_ref_end()
                for x in self.contextBamReads if (x.is_read2())]

    # Returns a list of all lengths for all BAM reads.
    def get_context_bam_read_lengths(self):
        return [x.get_bam_read_length() for x in self.contextBamReads]

    # Returns a list of BAM read sequences in the current context.
    def get_context_bam_read_seqs(self):
        return [x.get_bam_read_sequence() for x in self.contextBamReads]

    # Returns a list of qualities of all BAM reads.
    def get_context_bam_read_qualities(self):
        return [x.get_bam_read_qual() for x in self.contextBamReads]

    # Returns a list of Q-scores of all BAM reads.
    def get_context_bam_read_q_scores(self):
        return [x.get_bam_read_q_scores() for x in self.contextBamReads]

    # Returns a list of all BAM read MapQ values.
    def get_context_bam_read_map_qs(self):
        return [x.get_mapping_qual() for x in self.contextBamReads]

    # Returns whether a BAM read is in the context based on the provided
    # read identifier.
    def read_is_in_context(self, readId):
        return readId in self.get_context_bam_read_ids()

    # ===METHODS TO ADD/SET CONTEXT DATA=======================================
    # Adds the read id of a BAM read with an unmapped mate.
    def add_unmapped_mate_id(self, uReadId):
        self.unmappedReadMateIds.append(uReadId)

    # Sets the list of unmapped read mate ids.
    def set_unmapped_mate_ids(self, mateIds):
        self.unmappedReadMateIds = mateIds

    # Returns whether a BAM read in the context has an unmapped mate.
    def read_has_unmapped_mate(self, readId):
        return readId in self.unmappedReadMateIds

    # ===STATISTICS METHODS FOR A VARIANT CONTEXT==============================
    # Returns the average and median read length.
    def get_average_and_median_read_length(self):
        avgMedLen = []
        for contextread in self.contextBamReads:
            avgMedLen.append(contextread.get_bam_read_length())
        return ([statistics.mean(avgMedLen), statistics.median(avgMedLen)])

    # Returns the average and median read quality.
    def get_average_and_median_read_qual(self):
        avgMedQual = []
        for contextread in self.contextBamReads:
            avgMedQual.append(contextread.get_average_qscore())
        return ([statistics.mean(avgMedQual), statistics.median(avgMedQual)])

    # Returns the average and median read MapQ of this variant context.
    def get_average_and_median_read_map_q(self):
        avgMedMapQ = []
        for contextread in self.contextBamReads:
            avgMedMapQ.append(contextread.get_mapping_qual())
        return ([statistics.mean(avgMedMapQ), statistics.median(avgMedMapQ)])

    # ===SOME OTHER METHODS====================================================
    # Returns a string representation of the overlap context.
    def to_string(self):
        return (str(self.contextId) + "\t"
                + str(self.sampleId) + "\t"
                + str(self.contextChrom) + "\t"
                + str(self.contextOrigin) + "\t"
                + str(self.contextStart) + "\t"
                + str(self.contextEnd) + "\t"
                + str(len(self.contextBamReads)) + "\t"
                + ';'.join([x.get_bam_read_id() for x in self.contextBamReads]))

    # Returns a statistics string representation of the overlap context.
    def to_statistics_string(self):
        avgMedLens = self.get_average_and_median_read_length()
        avgMedQuals = self.get_average_and_median_read_qual()
        avgMedMapQ = self.get_average_and_median_read_map_q()
        return (f"{self.contextId}\t{avgMedLens[0]}\t{avgMedLens[1]}\t"
                f"{avgMedQuals[0]}\t{avgMedQuals[1]}\t{avgMedMapQ[0]}\t"
                f"{avgMedMapQ[1]}")

    # Compares the current OverlapContext to another OverlapContext and
    # returns the differences.
    def compare(self, otherOverlapContext):
        differences = {}
        if (self.contextId != otherOverlapContext.get_context_id()):
            differences[1] = [self.contextId,
                              otherOverlapContext.get_context_id()]
        if (self.sampleId != otherOverlapContext.get_sample_id()):
            differences[2] = [self.sampleId,
                              otherOverlapContext.get_sample_id()]
        if (self.contextChrom != otherOverlapContext.get_context_chrom()):
            differences[3] = [self.contextChrom,
                              otherOverlapContext.get_context_chrom()]
        if (self.contextOrigin != otherOverlapContext.get_context_origin()):
            differences[4] = [self.contextOrigin,
                              otherOverlapContext.get_context_origin()]
        if (self.contextStart != otherOverlapContext.get_context_start()):
            differences[5] = [self.contextStart,
                              otherOverlapContext.get_context_start()]
        if (self.contextEnd != otherOverlapContext.get_context_end()):
            differences[6] = [self.contextEnd,
                              otherOverlapContext.get_context_end()]
        if (self.get_context_bam_read_ids().sort() != otherOverlapContext.get_context_bam_read_ids().sort()):
            differences[7] = [self.contextBamReads,
                              otherOverlapContext.get_context_bam_reads()]
        return differences
