"""OverlapContext object class.

This module defines the OverlapContext object, which contains either the donor
or acceptor context information. These OverlapContext can be added to a
VariantContext object.
"""

import statistics as stats


class OverlapContext:
    """The OverlapContext saves an acceptor or donor context.

    The OverlapContext can be used to save an acceptor or donor context with
    reads that overlap with the variant directly.

    Attributes
    ----------
    context_id : str
        Identifier of the context
    sample_id : str
        Sample name/identifier the context was derived from
    context_chrom : str
        The chromosome name that the context is located on
    context_origin : int
        The variant genomic position that constructed the context
    context_start : int
        The leftmost genomic position of the context
    context_end : int
        The rightmost genomic position of the context
    context_reads : list of pysam.AlignedSegment
        List of aligned reads
    unmapped_read_mate_ids : list of str
        List of read IDs that have unmapped mates
    """

    def __init__(self, variantid, sampleid, ovconchrom, ovconorigin,
                 ovconstart, ovconend, bamreads):
        """Save the required data for an acceptor or donor context.

        Parameters
        ----------
        variantid : str
            Variant context identifier
        sampleid : str
            Sample name/identifier
        ovconchrom : str
            Chromosome name the context is located on
        ovconorigin : int
            Variant genomic position the context is constructed from
        ovconstart : int
            Leftmost genomic position of the context
        ovconend : int
            Rightmost genomic position of the context
        bamreads : list of pysam.AlignedSegment
            Reads and read mates overlapping with the context
        """
        self.context_id = variantid
        self.sample_id = sampleid
        self.context_chrom = ovconchrom
        self.context_origin = ovconorigin
        self.context_start = ovconstart
        self.context_end = ovconend
        self.context_reads = bamreads
        self.unmapped_read_mate_ids = []

    # ===METHODS TO GET DATA OF THE OVERLAP CONTEXT============================
    def get_context(self):
        """Return the essential data of the context.

        This data consists of the chromosome name the context is located on,
        the variant position the context is based on and the leftmost (start)
        and rightmost (end) genomic position of the context.

        Returns
        -------
        list
            List of essential context data
        """
        return [self.context_chrom, self.context_origin,
                self.context_start, self.context_end]

    def get_context_id(self):
        """Return the context identifier.

        Returns
        -------
        self.context_id : str
            The context identifier
        """
        return self.context_id

    def get_sample_id(self):
        """Return the sample identifier of the context.

        Returns
        -------
        self.sample_id : str
            The sample name/identifier
        """
        return self.sample_id

    def get_context_chrom(self):
        """Return the chromosome name that the context is located on.

        Returns
        -------
        self.context_chrom : str
            The context chromosome name
        """
        return self.context_chrom

    def get_context_origin(self):
        """Return the variant position that was used to construct the context.

        Returns
        -------
        self.context_origin : int
            Variant position the context is based on
        """
        return self.context_origin

    def get_context_start(self):
        """Return the left most genomic position of the context.

        Returns
        -------
        self.context_start : int
            Context starting position
        """
        return self.context_start

    def get_context_end(self):
        """Return the right most genomic position of the context.

        Returns
        -------
        self.context_end : int
            Context ending position
        """
        return self.context_end

    def get_context_bam_reads(self):
        """Return the list of reads and their mates overlapping with the context.

        Mates of reads overlapping with the context might not overlap with the
        context themselves. Each read is returned as a pysam AlignedSegment
        object.
        """
        return self.context_reads

    def get_unmapped_read_mate_ids(self):
        """Return a list with identifiers of reads that have an unmapped mate.

        Returns
        -------
        self.unmapped_read_mate_ids : list of str
            List with read IDs that have unmapped mates
        """
        return self.unmapped_read_mate_ids

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATION) FROM THE OVERLAP CONTEXT=============
    def get_context_length(self):
        """Calculate the distance between the variant context start and stop.

        Subtracts the left most genomic position of the context from the right
        most and returns the difference.

        Returns
        -------
        int
            Context length
        """
        return abs(self.context_end - self.context_start)

    def get_start_distance_from_origin(self):
        """Calculate the distance between the variant context start and variant position.

        Subtracts the left most position of the context from the variant
        position and returns the difference.

        Returns
        -------
        int
            Context start distance from variant position
        """
        return abs(self.context_origin - self.context_start)

    def get_end_distance_from_origin(self):
        """Calculate the distance between the variant position and variant context stop.

        Subtract the variant position from the rightmost position of the
        context and returns the difference.

        Returns
        -------
        int
            Context end distance form variant position
        """
        return abs(self.context_end - self.context_origin)

    # ===METHODS TO OBTAIN CONTEXT READ INFORMATION============================
    def get_number_of_context_reads(self):
        """Count the reads associated with the context.

        Returns
        -------
        int
            Number of context reads
        """
        if self.context_reads is None:
            return 0
        return len(self.context_reads)

    def get_context_read_ids(self):
        """Return a list of the identifiers of reads associated with the context.

        Returns
        -------
        list of str
            Identifiers of all reads
        """
        if self.context_reads is None:
            return [None]
        return list({x.query_name for x in self.context_reads})

    def get_context_read_starts(self):
        """Return a list of leftmost positions of all reads associated with the context.

        Returns
        -------
        list of int or None
            Leftmost positions of all reads
        """
        if self.context_reads is None:
            return [None]
        return [x.reference_start for x in self.context_reads]

    def get_context_read_left_positions(self):
        """Return a list of leftmost positions of all forward/R1 reads associated with the context.

        Returns
        -------
        list of int or None
            Leftmost genomic positions of all R1 reads
        """
        if self.context_reads is None:
            return [None]
        return [x.reference_start for x in self.context_reads if x.is_read1]

    def get_context_read_ends(self):
        """Return a list of rightmost positions of all reads associated with the context.

        Returns
        -------
        list of int or None
            Rightmost genomic positions of all reads
        """
        if self.context_reads is None:
            return [None]
        return [x.reference_end for x in self.context_reads]

    def get_context_read_right_positions(self):
        """Return a list of rightmost positions of all reverse/R2 reads associated with the context.

        Returns
        -------
        list of int
            Rightmost genomic positions of all reverse/R2 reads
        """
        if self.context_reads is None:
            return [None]
        return [x.reference_end for x in self.context_reads if x.is_read2]

    def get_context_read_lengths(self):
        """Return a list of lengths of all reads associated with the context.

        Returns
        -------
        list of int
            Lengths of all reads
        """
        return [x.reference_length for x in self.context_reads]

    def get_context_read_seqs(self):
        """Return a list of sequences of all reads associated with the context.

        Returns
        -------
        list of str
            Sequences of all reads
        """
        return [x.query_sequence for x in self.context_reads]

    def get_context_read_qualities(self):
        """Return a list of qualities lines of all reads associated with the context.

        Return
        ------
        list of str
            Quality line for all reads
        """
        read_qualities = []
        for cread in self.context_reads:
            read_qualities.append("".join([chr(x+33) for x in cread.query_qualities]))
        return read_qualities

    def get_context_read_q_scores(self):
        """Return a list of Q-Scores of all reads associated with the context.

        The Q-Scores for each read is an array of integers. The returned list
        is therefore a list with lists.

        Returns
        -------
        list
            Lists of Q-Scores for all reads
        """
        return [x.query_qualities for x in self.context_reads]

    def get_context_read_map_qs(self):
        """Return a list of MAPQ values of all reads.

        Returns
        -------
        list of int
            List with MAPQ values
        """
        return [x.mapping_quality for x in self.context_reads]

    def read_is_in_context(self, readid):
        """Check if a read specified by an ID is associated with the context.

        Returns
        -------
        bool
            True if the read in the context, False if not
        """
        if self.context_reads is None:
            return False
        return readid in self.get_context_read_ids()

    # ===METHODS TO ADD/SET CONTEXT DATA=======================================
    def add_unmapped_mate_id(self, ureadid):
        """Add the identifier of a read with an unmapped mate to the context unmapped mates list.

        Parameters
        ----------
        ureadid : str
            Identifier of the read with an unmapped mate
        """
        self.unmapped_read_mate_ids.append(ureadid)

    def set_unmapped_mate_ids(self, mateids):
        """Set the provided read ID list as the unmapped mate ID list associated with the context.

        Parameters
        ----------
        mateids : list of str
            read IDs with unmapped mates
        """
        self.unmapped_read_mate_ids = mateids

    def read_has_unmapped_mate(self, readid):
        """Check if a read specified by a read ID has an unmapped mate.

        Returns
        -------
        bool
            True if read has unmapped mate, False if not
        """
        return readid in self.unmapped_read_mate_ids

    # ===STATISTICS METHODS FOR A VARIANT CONTEXT==============================
    def get_average_and_median_read_length(self):
        """Calculate the mean and median read length of all reads associated with the context.

        Returns
        -------
        list of int
            Mean and median read length
        """
        if self.context_reads is None:
            return [None, None]
        avgmedlen = []
        for contextread in self.context_reads:
            if contextread.query_length is not None:
                avgmedlen.append(contextread.query_length)
        return [stats.mean(avgmedlen), stats.median(avgmedlen)]

    def get_average_and_median_read_qual(self):
        """Calculate the mean and median Q-Score of all reads associated with the context.

        Returns
        -------
        list of int
            Mean and median read Q-Score
        """
        if self.context_reads is None:
            return [None, None]
        avgmedqual = []
        for contextread in self.context_reads:
            avgmedqual.append(stats.mean(contextread.query_qualities))
        return [stats.mean(avgmedqual), stats.median(avgmedqual)]

    def get_average_and_median_read_map_q(self):
        """Calculate the mean and median MAPQ value of all reads associated with the context.

        Returns
        -------
        list of int
            Mean and median read MAPQ
        """
        if self.context_reads is None:
            return [None, None]
        avgmedmapq = []
        for contextread in self.context_reads:
            avgmedmapq.append(contextread.mapping_quality)
        return [stats.mean(avgmedmapq), stats.median(avgmedmapq)]

    # ===SOME OTHER METHODS====================================================
    def to_string(self):
        """Assemble and returns a string representation of the context.

        Returns
        -------
        str
            String representation of the context
        """
        if self.context_reads is None:
            bamids = None
            countbamreads = 0
        else:
            bamids = ";".join([x.query_name for x in self.context_reads])
            countbamreads = len(self.context_reads)
        return (str(self.context_id) + "\t"
                + str(self.sample_id) + "\t"
                + str(self.context_chrom) + "\t"
                + str(self.context_origin) + "\t"
                + str(self.context_start) + "\t"
                + str(self.context_end) + "\t"
                + str(countbamreads) + "\t"
                + str(bamids))

    def to_statistics_string(self):
        """Calculate some basic statistics of the context in a tab-separated string.

        Returns
        -------
        str
            Line with statistics to write to file
        """
        avgmedlens = self.get_average_and_median_read_length()
        avgmedquals = self.get_average_and_median_read_qual()
        avgmedmapq = self.get_average_and_median_read_map_q()
        return (f"{self.context_id}\t{avgmedlens[0]}\t{avgmedlens[1]}\t"
                f"{avgmedquals[0]}\t{avgmedquals[1]}\t{avgmedmapq[0]}\t"
                f"{avgmedmapq[1]}")

    def compare(self, other_overlap):
        """Compare the current context to a provided context and returns the differences.

        The context is compared to another context on each aspect. If they
        differ, the difference is saved in a dictionary. The key is a numeric
        value, the difference an array. The first entry is the value of the
        current context, the second entry is the value of the other context.

        Parameters
        ----------
        other_overlap : OverlapContext
            Acceptor/Donor context to compare this acceptor/donor context with

        Returns
        -------
        differences : dict
            Set of differences between the compared acceptor/donor contexts.
        """
        differences = {}
        if self.context_id != other_overlap.get_context_id():
            differences[1] = [self.context_id,
                              other_overlap.get_context_id()]
        if self.sample_id != other_overlap.get_sample_id():
            differences[2] = [self.sample_id,
                              other_overlap.get_sample_id()]
        if self.context_chrom != other_overlap.get_context_chrom():
            differences[3] = [self.context_chrom,
                              other_overlap.get_context_chrom()]
        if self.context_origin != other_overlap.get_context_origin():
            differences[4] = [self.context_origin,
                              other_overlap.get_context_origin()]
        if self.context_start != other_overlap.get_context_start():
            differences[5] = [self.context_start,
                              other_overlap.get_context_start()]
        if self.context_end != other_overlap.get_context_end():
            differences[6] = [self.context_end,
                              other_overlap.get_context_end()]
        if self.get_context_read_ids().sort() != other_overlap.get_context_read_ids().sort():
            differences[7] = [self.context_reads,
                              other_overlap.get_context_bam_reads()]
        return differences
