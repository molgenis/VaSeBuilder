import statistics


class OverlapContext:
    """The OverlapContext saves an acceptor or donor context.

    The OverlapContext can be used to save an acceptor or donor context with reads that overlap with the variant
    directly."""
    # Saves the data associated with the context overlapping with the variant.
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
    # Returns the context data (chrom, pos, start, end) in an array
    def get_context(self):
        """Returns the essential data of the context.

        This data consists of the chromosome name the context is located on, the variant position the context is based
        on and the leftmost (start) and rightmost (end) genomic position of the context."""
        return [self.context_chrom, self.context_origin, self.context_start, self.context_end]

    # Returns the context identifier.
    def get_context_id(self):
        """Returns the context identifier."""
        return self.context_id

    # Returns the sample identifier the context was based on.
    def get_sample_id(self):
        """Returns the identifier of the sample the context is """
        return self.sample_id

    # Returns the context chromosome name the context is located on.
    def get_context_chrom(self):
        """Returns the chromosome name of the context."""
        return self.context_chrom

    # Returns the origin position of the context.
    def get_context_origin(self):
        """Returns the the position of the variant that was used to construct the context."""
        return self.context_origin

    # Returns the context start position.
    def get_context_start(self):
        """Returns the left most genomic position of the context."""
        return self.context_start

    # Returns the end position of the context.
    def get_context_end(self):
        """Returns the right most genomic position of the context."""
        return self.context_end

    # Returns the bam reads associated with the context as a list of BamRead objects.
    def get_context_bam_reads(self):
        """Returns the list of reads and their mates overlapping with the context.

        Mates of reads overlapping with the context might not overlap with the context themselves. Each read is returned
        as a DonorBamRead object."""
        return self.context_bam_reads

    # Returns the list of BAM read ids that have an unmapped mate.
    def get_unmapped_read_mate_ids(self):
        """Returns a list with identifiers of reads that have an unmapped mate."""
        return self.unmapped_read_mate_ids

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATION) FROM THE OVERLAP CONTEXT=============
    # Returns the length of the context.
    def get_context_length(self):
        """Subtracts the left most genomic position of the context from the right most and returns the difference."""
        return abs(self.context_end - self.context_start)

    # Returns the distance of the context start from the context origin.
    def get_start_distance_from_origin(self):
        """Subtracts the left most position of the context from the variant position and returns the difference."""
        return abs(self.context_origin - self.context_start)

    # Returns the distance of the context end from the context origin.
    def get_end_distance_from_origin(self):
        """Subtracts the variant position from the rightmost position of the context and returns the difference."""
        return abs(self.context_end - self.context_origin)

    # ===METHODS TO OBTAIN CONTEXT READ INFORMATION============================
    # Returns the number of saved context reads.
    def get_number_of_context_reads(self):
        """Returns the number of reads associated with the context."""
        if self.context_bam_reads is None:
            return 0
        return len(self.context_bam_reads)

    # Returns a list of BAM read identifiers in the current context.
    def get_context_bam_read_ids(self):
        """Returns a list of the identifiers of reads associated with the context."""
        if self.context_bam_reads is None:
            return [None]
        return list(set([x.get_bam_read_id() for x in self.context_bam_reads]))

    # Returns a list of all left positions for all BAM reads.
    def get_context_bam_read_starts(self):
        """Returns a list of all let most position of all reads associated with the context."""
        if self.context_bam_reads is None:
            return [None]
        return [x.get_bam_read_ref_pos() for x in self.context_bam_reads]

    # Returns a list of all left positions for all R1 BAM reads.
    def get_context_bam_read_left_positions(self):
        """Returns a list with the left most positions for all forward/R1 reads associated with the context."""
        if self.context_bam_reads is None:
            return [None]
        return [x.get_bam_read_ref_pos()
                for x in self.context_bam_reads if (x.is_read1())]

    # Returns a list of BAM read ending positions for all BAM reads.
    def get_context_bam_read_ends(self):
        """Returns a list with the rightmost genomic positions of all reads associated with the context."""
        if self.context_bam_reads is None:
            return [None]
        return [x.get_bam_read_ref_end() for x in self.context_bam_reads]

    # Returns a list of all right positions for all R2 BAM reads.
    def get_context_bam_read_right_positions(self):
        """Returns a list with the rightmost genomic positions for all reverse/R2 reads associated with the context."""
        if self.context_bam_reads is None:
            return [None]
        return [x.get_bam_read_ref_end()
                for x in self.context_bam_reads if (x.is_read2())]

    # Returns a list of all lengths for all BAM reads.
    def get_context_bam_read_lengths(self):
        """Returns a list with the lengths of all reads associated with the context."""
        return [x.get_bam_read_length() for x in self.context_bam_reads]

    # Returns a list of BAM read sequences in the current context.
    def get_context_bam_read_seqs(self):
        """Returns a list with the sequences of all reads associated with the context."""
        return [x.get_bam_read_sequence() for x in self.context_bam_reads]

    # Returns a list of qualities of all BAM reads.
    def get_context_bam_read_qualities(self):
        return [x.get_bam_read_qual() for x in self.context_bam_reads]

    # Returns a list of Q-scores of all BAM reads.
    def get_context_bam_read_q_scores(self):
        """Returns a list with the Q-Scores of all reads associated with the context.

        The Q-Scores for each read is an array of integers. The returned list is therefore a list with lists."""
        return [x.get_bam_read_q_scores() for x in self.context_bam_reads]

    # Returns a list of all BAM read MapQ values.
    def get_context_bam_read_map_qs(self):
        return [x.get_mapping_qual() for x in self.context_bam_reads]

    # Returns whether a BAM read is in the context based on the provided read identifier.
    def read_is_in_context(self, readid):
        """Checks whether a read specified by an id is associated with the context."""
        if self.context_bam_reads is None:
            return False
        return readid in self.get_context_bam_read_ids()

    # ===METHODS TO ADD/SET CONTEXT DATA=======================================
    # Adds the read id of a BAM read with an unmapped mate.
    def add_unmapped_mate_id(self, ureadid):
        """Adds the identifier of a read with an unmapped mate to the context unmapped mates list."""
        self.unmapped_read_mate_ids.append(ureadid)

    # Sets the list of unmapped read mate ids.
    def set_unmapped_mate_ids(self, mateids):
        """Sets the provided read id list as the unmapped mate id list associated with the context."""
        self.unmapped_read_mate_ids = mateids

    # Returns whether a BAM read in the context has an unmapped mate.
    def read_has_unmapped_mate(self, readid):
        """Checks if a read specified by a read identifier has an unmapped mate and returns True or False."""
        return readid in self.unmapped_read_mate_ids

    # ===STATISTICS METHODS FOR A VARIANT CONTEXT==============================
    # Returns the average and median read length.
    def get_average_and_median_read_length(self):
        """Calculates and return the mean and median red length of all reads associated with the context."""
        if self.context_bam_reads is None:
            return [None, None]
        avgmedlen = []
        for contextread in self.context_bam_reads:
            if contextread.get_bam_read_length() is not None:
                avgmedlen.append(contextread.get_bam_read_length())
        return [statistics.mean(avgmedlen), statistics.median(avgmedlen)]

    # Returns the average and median read quality.
    def get_average_and_median_read_qual(self):
        """Calculates and returns the mean and median Q-Score of all reads associated with the context."""
        if self.context_bam_reads is None:
            return [None, None]        
        avgmedqual = []
        for contextread in self.context_bam_reads:
            avgmedqual.append(contextread.get_average_qscore())
        return [statistics.mean(avgmedqual), statistics.median(avgmedqual)]

    # Returns the average and median read MapQ of this variant context.
    def get_average_and_median_read_map_q(self):
        """Calculates the mean and median MAPQ value of all reads associated witht the context."""
        if self.context_bam_reads is None:
            return [None, None]
        avgmedmapq = []
        for contextread in self.context_bam_reads:
            avgmedmapq.append(contextread.get_mapping_qual())
        return [statistics.mean(avgmedmapq), statistics.median(avgmedmapq)]

    # ===SOME OTHER METHODS====================================================
    # Returns a string representation of the overlap context.
    def to_string(self):
        """Returns a String representation of the context."""
        if self.context_bam_reads is None:
            bamids = None
            countbamreads = 0
        else:
            bamids = ";".join([x.get_bam_read_id() for x in self.context_bam_reads])
            countbamreads = len(self.context_bam_reads)
        return (str(self.context_id) + "\t"
                + str(self.sample_id) + "\t"
                + str(self.context_chrom) + "\t"
                + str(self.context_origin) + "\t"
                + str(self.context_start) + "\t"
                + str(self.context_end) + "\t"
                + str(countbamreads) + "\t"
                + str(bamids))

    # Returns a statistics string representation of the overlap context.
    def to_statistics_string(self):
        """Calculates some basic statistics of the context in a tab separated String."""
        avgmedlens = self.get_average_and_median_read_length()
        avgmedquals = self.get_average_and_median_read_qual()
        avgmedmapq = self.get_average_and_median_read_map_q()
        return (f"{self.context_id}\t{avgmedlens[0]}\t{avgmedlens[1]}\t"
                f"{avgmedquals[0]}\t{avgmedquals[1]}\t{avgmedmapq[0]}\t"
                f"{avgmedmapq[1]}")

    # Compares the current OverlapContext to another OverlapContext and returns the differences.
    def compare(self, other_overlap_context):
        """Compares the current context to a provided context and returns the differences.

        The context is compared to another context on each aspect. If they differ, the difference is saved in a
        dictionary. The key is a numeric value, the difference an array. The first entry is the value of the current
        context, the second entry is the value of the other context."""
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
