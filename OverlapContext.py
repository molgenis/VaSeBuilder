import statistics
from DonorBamRead import DonorBamRead


class OverlapContext:
    # Saves the data associated with the context overlapping with the
    # variant.
    def __init__(self, variantid, sampleid, ovconchrom, ovconorigin,
                 ovconstart, ovconend, bamreads):
        self.context_id = variantid
        self.sample_id = sampleid
        self.context_chrom = ovconchrom
        self.context_origin = ovconorigin
        self.context_start = ovconstart
        self.context_end = ovconend
        self.context_bam_reads = bamreads
        self.unmapped_read_mate_ids = []

    # ===METHODS TO GET DATA OF THE OVERLAP CONTEXT============================
    # Returns the context identifier.
    def get_context_id(self):
        return self.context_id

    # Returns the sample identifier the context was based on.
    def get_sample_id(self):
        return self.sample_id

    # Returns the context chromosome name the context is located on.
    def get_context_chrom(self):
        return self.context_chrom

    # Returns the origin position of the context.
    def get_context_origin(self):
        return self.context_origin

    # Returns the context start position.
    def get_context_start(self):
        return self.context_start

    # Returns the end position of the context.
    def get_context_end(self):
        return self.context_end

    # Returns the bam reads associated with the context as a list of
    # BamRead objects.
    def get_context_bam_reads(self):
        return self.context_bam_reads

    # Returns the list of BAM read ids that have an unmapped mate.
    def get_unmapped_read_mate_ids(self):
        return self.unmapped_read_mate_ids

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATION) FROM THE=============
    # ===OVERLAP CONTEXT===
    # Returns the lengt of the context.
    def get_context_length(self):
        return abs(self.context_end - self.context_start)

    # Returns the distance of the context start from the context origin.
    def get_start_distance_from_origin(self):
        return abs(self.context_origin - self.context_start)

    # Returns the distance of the context end from the context origin.
    def get_end_distance_from_origin(self):
        return abs(self.context_end - self.context_origin)

    # ===METHODS TO OBTAIN CONTEXT READ INFORMATION============================
    # Returns the number of saved context reads.
    def get_number_of_context_reads(self):
        return len(self.context_bam_reads)

    # Returns a list of BAM read identifiers in the current context.
    def get_context_bam_read_ids(self):
        return [x.get_bam_read_id() for x in self.context_bam_reads]

    # Returns a list of all left positions for all BAM reads.
    def get_context_bam_read_starts(self):
        return [x.get_bam_read_ref_pos() for x in self.context_bam_reads]

    # Returns a list of all left positions for all R1 BAM reads.
    def get_context_bam_read_left_positions(self):
        return [x.get_bam_read_ref_pos()
                for x in self.context_bam_reads if (x.is_read1())]

    # Returns a list of BAM read ending positions for all BAM reads.
    def get_context_bam_read_ends(self):
        return [x.get_bam_read_ref_end() for x in self.context_bam_reads]

    # Returns a list of all right positions for all R2 BAM reads.
    def get_context_bam_read_right_positions(self):
        return [x.get_bam_read_ref_end()
                for x in self.context_bam_reads if (x.is_read2())]

    # Returns a list of all lengths for all BAM reads.
    def get_context_bam_read_lengths(self):
        return [x.get_bam_read_length() for x in self.context_bam_reads]

    # Returns a list of BAM read sequences in the current context.
    def get_context_bam_read_seqs(self):
        return [x.get_bam_read_sequence() for x in self.context_bam_reads]

    # Returns a list of qualities of all BAM reads.
    def get_context_bam_read_qualities(self):
        return [x.get_bam_read_qual() for x in self.context_bam_reads]

    # Returns a list of Q-scores of all BAM reads.
    def get_context_bam_read_q_scores(self):
        return [x.get_bam_read_q_scores() for x in self.context_bam_reads]

    # Returns a list of all BAM read MapQ values.
    def get_context_bam_read_map_qs(self):
        return [x.get_mapping_qual() for x in self.context_bam_reads]

    # Returns whether a BAM read is in the context based on the provided
    # read identifier.
    def read_is_in_context(self, readid):
        return readid in self.get_context_bam_read_ids()

    # ===METHODS TO ADD/SET CONTEXT DATA=======================================
    # Adds the read id of a BAM read with an unmapped mate.
    def add_unmapped_mate_id(self, ureadid):
        self.unmapped_read_mate_ids.append(ureadid)

    # Sets the list of unmapped read mate ids.
    def set_unmapped_mate_ids(self, mateids):
        self.unmapped_read_mate_ids = mateids

    # Returns whether a BAM read in the context has an unmapped mate.
    def read_has_unmapped_mate(self, readid):
        return readid in self.unmapped_read_mate_ids

    # ===STATISTICS METHODS FOR A VARIANT CONTEXT==============================
    # Returns the average and median read length.
    def get_average_and_median_read_length(self):
        avgmedlen = []
        for contextread in self.context_bam_reads:
            avgmedlen.append(contextread.get_bam_read_length())
        return [statistics.mean(avgmedlen), statistics.median(avgmedlen)]

    # Returns the average and median read quality.
    def get_average_and_median_read_qual(self):
        avgmedqual = []
        for contextread in self.context_bam_reads:
            avgmedqual.append(contextread.get_average_qscore())
        return [statistics.mean(avgmedqual), statistics.median(avgmedqual)]

    # Returns the average and median read MapQ of this variant context.
    def get_average_and_median_read_map_q(self):
        avgmedmapq = []
        for contextread in self.context_bam_reads:
            avgmedmapq.append(contextread.get_mapping_qual())
        return [statistics.mean(avgmedmapq), statistics.median(avgmedmapq)]

    # ===SOME OTHER METHODS====================================================
    # Returns a string representation of the overlap context.
    def to_string(self):
        return (str(self.context_id) + "\t"
                + str(self.sample_id) + "\t"
                + str(self.context_chrom) + "\t"
                + str(self.context_origin) + "\t"
                + str(self.context_start) + "\t"
                + str(self.context_end) + "\t"
                + str(len(self.context_bam_reads)) + "\t"
                + ";".join([x.get_bam_read_id() for x in self.context_bam_reads]))

    # Returns a statistics string representation of the overlap context.
    def to_statistics_string(self):
        avgmedlens = self.get_average_and_median_read_length()
        avgmedquals = self.get_average_and_median_read_qual()
        avgmedmapq = self.get_average_and_median_read_map_q()
        return (f"{self.context_id}\t{avgmedlens[0]}\t{avgmedlens[1]}\t"
                f"{avgmedquals[0]}\t{avgmedquals[1]}\t{avgmedmapq[0]}\t"
                f"{avgmedmapq[1]}")

    # Compares the current OverlapContext to another OverlapContext and
    # returns the differences.
    def compare(self, other_overlap_context):
        differences = {}
        if self.context_id != other_overlap_context.get_context_id():
            differences[1] = [self.context_id,
                              other_overlap_context.get_context_id()]
        if self.sample_id != other_overlap_context.get_sample_id():
            differences[2] = [self.sample_id,
                              other_overlap_context.get_sample_id()]
        if self.context_chrom != other_overlap_context.get_context_chrom():
            differences[3] = [self.context_chrom,
                              other_overlap_context.get_context_chrom()]
        if self.context_origin != other_overlap_context.get_context_origin():
            differences[4] = [self.context_origin,
                              other_overlap_context.get_context_origin()]
        if self.context_start != other_overlap_context.get_context_start():
            differences[5] = [self.context_start,
                              other_overlap_context.get_context_start()]
        if self.context_end != other_overlap_context.get_context_end():
            differences[6] = [self.context_end,
                              other_overlap_context.get_context_end()]
        if self.get_context_bam_read_ids().sort() != other_overlap_context.get_context_bam_read_ids().sort():
            differences[7] = [self.context_bam_reads,
                              other_overlap_context.get_context_bam_reads()]
        return differences
