"""Variant Context object class.

This module creates the VariantContext object, which stores all information
related to each individual variant context. It also defines numerous methods
related to its functionality.
"""

import statistics
from overlap_context import OverlapContext


class VariantContext:
    """Save a variant context and allows data to be obtained and calculated.

    Attributes
    ----------
    context_id : str
        Identifier of the context
    sample_id : str
        The ID of the donor sample used to construct the variant context
    variant_context_chrom : str
        Chromosome name the context is lcoated on
    variant_context_origin : int
        Variant position the context is constructed from
    variant_context_start : int
        Leftmost genomic position of the variant context
    variant_context_end : int
        Rightmost genomic position of the variant context
    variant_context_areads : list of pysam.AlignedSegment
        Acceptor reads associated with the variant context
    variant_context_dreads : list of pysam.AlignedSegment
        Donor reads associated with the variant context
    variant_acceptor_context : OverlapContext
        Acceptor context used to construct the variant context
    variant_donor_context : OverlapContext
        Donor context used to construct the variant context
    unmapped_donor_mate_ids : list of str
        IDs of reads with an unmapped mate
    """

    def __init__(self, varconid, sampleid,
                 varconchrom, varconorigin,
                 varconstart, varconend,
                 acceptorreads, donorreads,
                 acceptor_context=None, donor_context=None,
                 variants=None):
        """Save the provided variant context data.

        Acceptor and donor contexts are optional when creating the variant
        context.

        Parameters
        ----------
        varconid : str
            Variant context identifier
        sampleid : str
            Sample name/identifier
        varconchrom : str
            Chromosome name the variant context is located on
        varconorigin:
            Variant position the context is constructed from
        varconstart : int
            Leftmost genomic position of the variant context
        varconend : int
            Rightmost genomic position of the variant context
        acceptorreads : list of pysam.AlignedSegment
            Variant context acceptor reads
        donorreads : list of pysam.AlignedSegment
            Variant context donor reads
        acceptor_context : OverlapContext
            Acceptor context associated with the variant context
        donor_context : OverlapContext
            Donor context associated with the variant context
        """
        self.context_id = varconid
        self.sample_id = sampleid
        self.variant_context_chrom = varconchrom
        self.variant_context_origin = int(varconorigin)
        self.variant_context_start = int(varconstart)
        self.variant_context_end = int(varconend)
        # Acceptor reads that overlap with the entire variant context.
        self.variant_context_areads = acceptorreads
        # Donor reads that overlap with the entire variant context.
        self.variant_context_dreads = donorreads
        # Context from acceptor reads overlapping with the variant itself.
        self.variant_acceptor_context = acceptor_context
        # Context from donor reads overlapping with the variant itself.
        self.variant_donor_context = donor_context
        self.unmapped_acceptor_mate_ids = []
        self.unmapped_donor_mate_ids = []
        # self.priority_label = None
        # self.priority_level = None
        self.priorities = []
        self.variants = variants
        if self.variants is None:
            self.variants = []

    # ===METHODS TO OBTAIN DATA FROM THE VARIANT CONTEXT DATA==================
    def get_context(self):
        """Return the essential data of the variant context.

        Returns
        -------
        list of str and int
            Variant context window
        """
        return [self.variant_context_chrom, self.variant_context_origin,
                self.variant_context_start, self.variant_context_end]

    def get_variant_context_id(self):
        """Return the variant context identifier.

        Returns
        -------
        self.context_id : str
            Variant context identifier
        """
        return self.context_id

    def get_variant_context_sample(self):
        """Return the name/identifier of the sample.

        Returns
        -------
        self.sample_id : str
            Sample name/identifier
        """
        return self.sample_id

    def get_variant_context_chrom(self):
        """Return the chromosome name the variant context is located on.

        Returns
        -------
        self.variant_context_chrom : str
            Chromosome name of the variant context
        """
        return self.variant_context_chrom

    def get_variant_context_origin(self):
        """Return the position of the variant used to construct the variant context.

        Returns
        -------
        self.variant_context_origin : int
            Variant genomic position
        """
        return self.variant_context_origin

    def get_variant_context_start(self):
        """Return the leftmost genomic position of the variant context.

        Returns
        -------
        self.variant_context_start
            Variant context starting position
        """
        return self.variant_context_start

    def get_variant_context_end(self):
        """Return the variant context end position.

        Returns
        -------
        self.variant_context_end : int
            Variant context rightmost genomic position
        """
        return self.variant_context_end

    def get_acceptor_reads(self):
        """Return the variant context acceptor reads.

        Returns
        -------
        self.variant_context_areads : list of pysam.AlignedSegment
            Variant context acceptor reads
        """
        return self.variant_context_areads

    def get_donor_reads(self):
        """Return the variant context donor reads.

        Returns
        -------
        self.variant_context_dreads : list of pysam.AlignedSegment
            Variant context acceptor reads
        """
        return self.variant_context_dreads

    def set_donor_reads(self, donor_reads):
        """Set a list of reads as the variant context donor reads.

        Parameters
        ----------
        donor_reads: list of pysam.AlignedSegment
            Donor reads to set
        """
        self.variant_context_dreads = donor_reads

    def get_donor_read_strings(self):
        """Return all donor reads as tuples of strings.

        Tuples each contain 4 strings: Read name, pair number, sequence,
        and quality line.

        Returns
        -------
        list(set(donorreads))
            List of donor read string tuples.
        """
        donorreads = []
        for dbr in self.variant_context_dreads:
            readpn = "2"
            if dbr.is_read1:
                readpn = "1"
            donorreads.append(
                (dbr.query_name,
                 readpn,
                 dbr.query_sequence,
                 "".join([chr(x+33) for x in dbr.query_qualities]))
                )
        return list(set(donorreads))

    def get_acceptor_context(self):
        """Return the acceptor context used to construct the variant context.

        Returns
        -------
        self.variant_acceptor_context : OverlapContext
            Acceptor context associated with the variant context
        """
        return self.variant_acceptor_context

    def get_donor_context(self):
        """Return the donor context used to construct the variant context.

        Returns
        -------
        self.variant_donor_context : OverlapContext
            Donor context associated with the variant context
        """
        return self.variant_donor_context

    def get_unmapped_acceptor_mate_ids(self):
        """Return the IDs of variant context acceptor reads with unmapped mates.

        Returns
        -------
        self.unmapped_acceptor_mate_ids : list of str
            Variant context acceptor read IDs with unmapped mates
        """
        return self.unmapped_acceptor_mate_ids

    # Returns a list of donor reads overlapping with the variant context.
    def get_unmapped_donor_mate_ids(self):
        """Return the IDs of variant context donor reads with unmapped mates.

        Returns
        -------
        self.unmapped_donor_mate_ids : list of str
            Variant contexct donor read IDs with an unmapped mate
        """
        return self.unmapped_donor_mate_ids

# =============================================================================
#     def get_priority_label(self):
#         """Return the priority label of the variant context.
#
#         Returns
#         -------
#         self.priority_label : str
#             The priority label associated with the variant context
#         """
#         return self.priority_label
#
#     def get_priority_level(self):
#         """Return the priority level of the variant context.
#
#         Returns
#         -------
#         self.priority_level : int
#             Priority level of the variant context
#         """
#         return self.priority_level
#
#     def set_priority_label(self, prlabel):
#         """Set the priority label for the variant context.
#
#         Parameters
#         ----------
#         prlabel : str
#             Priority label to set for the variant context
#         """
#         self.priority_label = prlabel
#
#     def set_priority_level(self, prlevel):
#         """Set the priority level for the variant context.
#
#         Parameters
#         ----------
#         prlevel : int
#             Priority level to set for the variant context
#         """
#         self.priority_level = prlevel
# =============================================================================

    # ===METHODS TO GET CALCULATED DATA OF THE VARIANT CONTEXT=================
    def get_variant_context_length(self):
        """Calculate the length of the variant context.

        The length if determined by subtracting the leftmost genomic (start)
        position from the rightmost genomic (end) position.

        Returns
        -------
        int
            Variant context length
        """
        return abs(self.variant_context_end - self.variant_context_start)

    def get_start_distance_from_origin(self):
        """Calculate the distance from variant context start to variant locus.

        Returns
        -------
        int
            Distance between variant context start and variant start positions
        """
        return abs(self.variant_context_origin - self.variant_context_start)

    def get_end_distance_from_origin(self):
        """Calculate the distance from variant locus to variant context end.

        Returns
        -------
        int
            Distance between variant start and variant context end positions
        """
        return abs(self.variant_context_end - self.variant_context_origin)

    # ===METHODS TO OBTAIN VARIANT CONTEXT ACCEPTOR READ DATA==================
    def get_number_of_acceptor_reads(self):
        """Return the number of variant context acceptor reads.

        Returns
        -------
        int
            Number of variant context acceptor reads
        """
        if self.variant_context_areads is None:
            return 0
        return len(self.variant_context_areads)

    def get_acceptor_read_ids(self):
        """Return the variant context acceptor read IDs.

        Returns
        -------
        list of str or None
            Variant context acceptor read IDs, None if there are none
        """
        if self.variant_context_areads is None:
            return [None]
        return list({x.query_name for x in self.variant_context_areads})

    def get_acceptor_read_starts(self):
        """Return the variant context acceptor read starting positions.

        Returns
        -------
        list of int or None
            Variant context acceptor read leftmost genomic positions,
            None if there are no acceptor reads
        """
        if self.variant_context_areads is None:
            return [None]
        return [x.reference_start for x in self.variant_context_areads]

    def get_acceptor_read_left_positions(self):
        """Return the leftmost genomic positions of all variant context R1 acceptor reads.

        Returns
        -------
        list of int or None
            Variant context R1 acceptor read leftmost genomic positions,
            None if there are no acceptor reads
        """
        if self.variant_context_areads is None:
            return [None]
        return [x.reference_start for x in self.variant_context_areads if x.is_read1]

    def get_acceptor_read_ends(self):
        """Return the variant context acceptor read rightmost positions.

        Returns
        -------
        list of int or None
            Variant context R2 acceptor read rightmost genomic positions,
            None if there are no acceptor reads
        """
        if self.variant_context_areads is None:
            return [None]
        return [x.reference_end for x in self.variant_context_areads]

    def get_acceptor_read_right_positions(self):
        """Return the rightmost genomic positions for all variant context R2 acceptor reads.

        Returns
        -------
        list of int or None
            Variant context
        """
        if self.variant_context_areads is None:
            return [None]
        return [x.reference_end for x in self.variant_context_areads
                if x.is_read2]

    # ===METHODS TO OBTAIN VARIANT CONTEXT DONOR READ DATA=====================
    def get_number_of_donor_reads(self):
        """Return the number of variant context donor reads.

        Returns
        -------
        int
            Number of variant context donor reads
        """
        return len(self.variant_context_dreads)

    def get_donor_read_ids(self):
        """Return the IDs of donor reads overlapping with the variant context.

        Returns
        -------
        list of str
            Variant context donor read IDs
        """
        return list({x.query_name for x in self.variant_context_dreads})

    def get_donor_read_starts(self):
        """Return the list of variant context donor read starting positions.

        Returns
        -------
        list of int
            Variant context donor read leftmost genomic positions
        """
        return [x.reference_start for x in self.variant_context_dreads]

    def get_donor_read_left_positions(self):
        """Return the list of variant context R1 donor read leftmost positions.

        Returns
        -------
        list of int
            Variant context R1 donor read leftmost genomic positions
        """
        return [x.reference_start for x in self.variant_context_dreads
                if x.is_read1]

    # Returns a list of all donor read ending positions.
    def get_donor_read_ends(self):
        """Return variant context donor read rightmost positions.

        Returns
        -------
        list of int
            Variant context donor read rightmost genomic positions
        """
        return [x.reference_end for x in self.variant_context_dreads]

    def get_donor_read_right_positions(self):
        """Return the list of variant context R2 donor reads.

        Returns
        -------
        list of int
            Variant context R2 donor reads rightmost genomic positions
        """
        return [x.reference_end for x in self.variant_context_dreads
                if x.is_read2]

    # ===METHODS TO ADD DATA TO THE VARIANT CONTEXT============================
    def set_acceptor_context(self, acceptor_context):
        """Set the acceptor context from an existing context.

        Parameters
        ----------
        acceptor_context : OverlapContext
            Acceptor context to add to the variant context
        """
        self.variant_acceptor_context = acceptor_context

    def set_donor_context(self, donor_context):
        """Set the donor context of the variant context with the one provided.

        Parameters
        ----------
        donor_context : OverlapContext
            Donor context to add to the variant context
        """
        self.variant_donor_context = donor_context

    def add_acceptor_context(self, contextid, sampleid,
                             contextchrom, contextorigin,
                             contextstart, contextend,
                             acceptorreads):
        """Construct and set the acceptor context from the provided parameters.

        Parameters
        ----------
        contextid : str
            Acceptor context identifier
        sampleid : str
            Sample name/identifier
        contextchrom : str
            Chromosome name the context is located on
        contextorigin : int
            Variant position the context is constructed from
        contextstart : int
            Leftmost genomic position of the acceptor context
        contextend : int
            Rightmost genomic position of the acceptor context
        acceptorreads : list of pysam.AlignedSegment
            Acceptor context reads
        """
        self.variant_acceptor_context = OverlapContext(
            contextid, sampleid,
            contextchrom, contextorigin,
            contextstart, contextend,
            acceptorreads
            )

    def add_donor_context(self, contextid, sampleid,
                          contextchrom, contextorigin,
                          contextstart, contextend,
                          donorreads):
        """Construct and set the donor context from the provided parameters.

        Parameters
        ----------
        contextid : str
            Donor context identifier
        sampleid : str
            Sample name/identifier
        contextchrom : str
            Chromosome name the context is located on
        contextorigin : int
            Variant genomic position the context is constructed from
        contextstart : int
            Leftmost genomic position of the donor context
        contextend : int
            Rightmost genomic position of the donor context
        donorreads : list of pysam.AlignedSegment
            Donor context reads
        """
        self.variant_donor_context = OverlapContext(
            contextid, sampleid,
            contextchrom, contextorigin,
            contextstart, contextend,
            donorreads
            )

    # ===METHODS TO OBTAIN VARIANT CONTEXT UNMAPPED MATE READ DATA=============
    def get_unmapped_acceptor_read_ids(self):
        """Return the read IDs of variant context acceptor reads with unmapped mates.

        Returns
        -------
        self.unmapped_acceptor_mate_ids : list of str
            Variant context acceptor read IDs with an unmapped mate
        """
        return self.unmapped_acceptor_mate_ids

    def get_unmapped_donor_read_ids(self):
        """Return the read IDs of variant context donor reads with unmapped mates.

        Returns
        -------
        self.unmapped_donor_mate_ids : list of str
            Variant context donor read IDs with an unmapped mate
        :return:
        """
        return self.unmapped_donor_mate_ids

    def add_unmapped_acceptor_mate_id(self, mateid):
        """Add a variant context appector mate ID.

        Parameters
        ----------
        mateid : str
            Variant context acceptor read ID with an unmapped mate
        """
        self.unmapped_acceptor_mate_ids.append(mateid)

    def add_unmapped_donor_mate_id(self, mateid):
        """Add a variant context donor mate ID.

        Parameters
        ----------
        mateid : str
            Variant context donor read ID with an unmapped mate
        """
        self.unmapped_donor_mate_ids.append(mateid)

    def set_unmapped_acceptor_mate_ids(self, mateids):
        """Set the variant context unmapped acceptor mate ids.

        Parameters
        ----------
        mateids : list of str
            Variant context acceptor read IDs with unmapped mate
        """
        self.unmapped_acceptor_mate_ids = mateids

    def set_unmapped_donor_mate_ids(self, mateids):
        """Set the variant context unmapped donor mate ids.

        Parameters
        ----------
        mateids : list of str
            Variant context donor read IDs with unmapped mates
        """
        self.unmapped_donor_mate_ids = mateids

    def acceptor_read_has_unmapped_mate(self, readid):
        """Check if a specified variant context acceptor read has an unmapped mate.

        Parameters
        ----------
        readid : str
            Acceptor read ID

        Returns
        -------
        bool
            True if acceptor read has unmapped mate, False if not
        """
        return readid in self.unmapped_acceptor_mate_ids

    def donor_read_has_unmapped_mate(self, readid):
        """Check if a specified variant context donor read has an unmapped mate.

        Parameters
        ----------
        readid : str
            Donor read ID

        Returns
        -------
        bool
            True if donor read has unmapped mate, False if not
        """
        return readid in self.unmapped_donor_mate_ids

    def get_number_of_unmapped_acceptor_mates(self):
        """Return the number of variant context acceptor reads with unmapped mates.

        Returns
        -------
        int
            Number of variant context acceptor reads with an unmapped mate.
        """
        return len(self.unmapped_acceptor_mate_ids)

    def get_number_of_unmapped_donor_mates(self):
        """Return the number of variant context donor reads with unmapped mates.

        Returns
        -------
        int
            Number of variant context donor reads with an unmapped mate.
        """
        return len(self.unmapped_donor_mate_ids)

    # ===METHODS TO ADD UNMAPPED MATES TO THE ACCEPTOR AND DONOR CONTEXT=======
    def set_acceptor_context_unmapped_mates(self, mateids):
        """Set the unmapped mate ids for the acceptor context.

        Returns
        -------
        mateids : list of str
            Acceptor context read IDs with an unmapped mate
        """
        self.variant_acceptor_context.set_unmapped_mate_ids(mateids)

    def add_acceptor_context_unmapped_mate(self, ureadid):
        """Add an unmapped read id to the acceptor context.

        Parameters
        ----------
        ureadid : str
            Acceptor context read ID with an unmapped mate
        """
        self.variant_acceptor_context.add_unmapped_mate_id(ureadid)

    def set_donor_context_unmapped_mates(self, mateids):
        """Set the unmapped mate ids for the donor context.

        Parameters
        ----------
        mateids : list of str
            Donor context read IDs with an unmapped mate
        """
        self.variant_donor_context.set_unmapped_mate_ids(mateids)

    def add_donor_context_unmapped_mate(self, ureadid):
        """Add an unmapped read id to the donor context.

        Parameters
        ----------
        ureadid : str
            Donor context read ID with an unmapped mate
        """
        self.variant_donor_context.add_unmapped_mate_id(ureadid)

    # ===METHODS TO OBTAIN STATISTICS OF THE VARIANT CONTEXT===================
    def get_average_and_median_acceptor_read_qual(self):
        """Calculate the mean and median variant context acceptor read quality.

        Returns
        -------
        qual_stats : list of float or list of None
            Average and median read quality
        """
        qual_stats = self.get_average_and_median_read_qual(self.variant_context_areads)
        return qual_stats

    def get_average_and_median_donor_read_qual(self):
        """Calculate the mean and median variant context donor read quality.

        Returns
        -------
        qual_stats : list of float or list of None
            Average and median read quality
        """
        qual_stats = self.get_average_and_median_read_qual(self.variant_context_dreads)
        return qual_stats

    @staticmethod
    def get_average_and_median_read_qual(contextreads):
        """Calculate the mean and median read quality score.

        Parameters
        ----------
        contextreads : list or reads
            Reads to calculate mean and median quality score of

        Returns
        -------
        list of float or list of None
            Average and median read quality
        """
        if contextreads is not None:
            avgmedqual = []
            for contextread in contextreads:
                avgmedqual.append(statistics.mean(list(contextread.query_qualities)))
            return ([statistics.mean(avgmedqual),
                     statistics.median(avgmedqual)])
        return [None, None]

    def get_average_and_median_acceptor_read_mapq(self):
        """Calculate the mean and median variant context acceptor read MAPQ values.

        Returns
        -------
        list of int
            Mean and median MAPQ, None if there are no acceptor reads
        """
        return self.get_average_and_median_read_mapq(self.variant_context_areads)

    def get_average_and_median_donor_read_mapq(self):
        """Calculate the mean and median variant context donor read MAPQ values.

        Returns
        -------
        list of int
            Mean and median variant context donor read MAPQ, None if there are no donor reads
        """
        return self.get_average_and_median_read_mapq(self.variant_context_dreads)

    @staticmethod
    def get_average_and_median_read_mapq(contextreads):
        """Calculate the mean and median MAPQ value of provided reads.

        Parameters
        ----------
        contextreads : list of pysam.AlignedSegment
            Reads to calculate mean and median MAPQ of

        Returns
        -------
        list of int
            Mean and median read MAPQ, None if no reads provided
        """
        if contextreads is not None:
            avgmedmapq = []
            for contextread in contextreads:
                avgmedmapq.append(contextread.mapping_quality)
            return ([statistics.mean(avgmedmapq),
                     statistics.median(avgmedmapq)])
        return [None, None]

    def get_average_and_median_acceptor_read_length(self):
        """Calculate the mean and median variant context acceptor read length.

        Returns
        -------
        list of int
            Mean and median variant context acceptor read length, None if
            there are no acceptor reads
        """
        return self.get_average_and_median_read_length(self.variant_context_areads)

    def get_average_and_median_donor_read_length(self):
        """Calculate the mean and median variant context donor read length.

        Returns
        -------
        list of int
            Mean and median variant context donor read length, None if there are no donor reads
        """
        return self.get_average_and_median_read_length(self.variant_context_dreads)

    @staticmethod
    def get_average_and_median_read_length(contextreads):
        """Calculate the mean and median read length of a specified list of reads.

        Parameters
        ----------
        contextreads : list of pysam.AlignedSegment
            Reads to to calculate mean and median length of

        Returns
        -------
        list of int
            Mean and median read length, None if no reads are supplied
        """
        if contextreads is not None:
            avgmedlen = []
            for contextread in contextreads:
                if contextread.reference_length is not None:
                    avgmedlen.append(contextread.reference_length)
            return [statistics.mean(avgmedlen), statistics.median(avgmedlen)]
        return [None, None]

    # ===METHODS TO OBTAIN ACCEPTOR CONTEXT DATA===============================
    def has_acceptor_context(self):
        """Return whether the variant context has an acceptor context.

        Returns
        -------
        bool
            True if variant context has an acceptor context, False if not
        """
        return self.variant_acceptor_context is not None

    def get_acceptor_context_id(self):
        """Return the acceptor context identifier.

        Returns
        -------
        str
            Acceptor context identifier
        """
        return self.variant_acceptor_context.get_context_id()

    def get_acceptor_context_sample_id(self):
        """Return the acceptor context sample id.

        Returns
        -------
        str
            Sample name/identifier of the acceptor context
        """
        return self.variant_acceptor_context.get_sample_id()

    def get_acceptor_context_chrom(self):
        """Return the chromosome name of the acceptor context.

        Returns
        -------
        str
            Chromosome name the acceptor context is located on
        """
        return self.variant_acceptor_context.get_context_chrom()

    def get_acceptor_context_origin(self):
        """Return the origin position of the acceptor context.

        Returns
        -------
        int
            Variant genomic position the context is based on
        """
        return self.variant_acceptor_context.get_context_origin()

    def get_acceptor_context_start(self):
        """Return the starting position of the acceptor context.

        Returns
        -------
        int
            Acceptor context leftmost genomic position
        """
        return self.variant_acceptor_context.get_context_start()

    def get_acceptor_context_end(self):
        """Return the ending position of the acceptor context.

        Returns
        -------
        int
            Acceptor context rightmost genomic position
        """
        return self.variant_acceptor_context.get_context_end()

    def get_acceptor_context_length(self):
        """Return the length of the acceptor context.

        Returns
        -------
        int
            Acceptor context length
        """
        return self.variant_acceptor_context.get_context_length()

    def get_acceptor_context_reads(self):
        """Return the acceptor context reads.

        Returns
        -------
        list of pysam.AlignedSegment
            Acceptor context reads
        """
        return self.variant_acceptor_context.get_context_bam_reads()

    def get_acceptor_context_read_ids(self):
        """Return the acceptor context read IDs.

        Returns
        -------
        list of str
            Acceptor context read IDs
        """
        return self.variant_acceptor_context.get_context_read_ids()

    def get_acceptor_context_read_starts(self):
        """Return the leftmost genomic positions of the acceptor context reads.

        Returns
        -------
        list of int
            List of read start positions
        """
        return self.variant_acceptor_context.get_context_read_starts()

    def get_acceptor_context_read_left_positions(self):
        """Return the leftmost genomic positions of all R1 acceptor context reads.

        Returns
        -------
        list of int
            Acceptor context leftmost genomic R1 read positions
        """
        return self.variant_acceptor_context.get_context_read_left_positions()

    def get_acceptor_context_read_ends(self):
        """Return the rightmost genomic positions of all acceptor context reads.

        Returns
        -------
        list of int
            Acceptor context rightmost genomic read positions.
        """
        return self.variant_acceptor_context.get_context_read_ends()

    def get_acceptor_context_read_right_positions(self):
        """Return a list of all acceptor context R2 BAM read end positions.

        Returns
        -------
        list of int
            Acceptor context rightmost genomic R2 read positions
        """
        return self.variant_acceptor_context.get_context_read_right_positions()

    def get_acceptor_context_read_lengths(self):
        """Return the lengths of the acceptor context reads.

        Returns
        -------
        list of int
            Acceptor context read lengths
        """
        return self.variant_acceptor_context.get_context_read_lengths()

    def get_acceptor_context_unmapped_mate_ids(self):
        """Return the acceptor context read IDs with unmapped mates.

        Returns
        -------
        list of str
            Acceptor context read IDs with unmapped mates.
        """
        return self.variant_acceptor_context.get_unmapped_read_mate_ids()

    # ===METHODS TO OBTAIN DONOR CONTEXT DATA==================================
    def has_donor_context(self):
        """Check whether the variant context has a donor context saved.

        Returns
        -------
        bool
            True if the variant context has a donor context, False if not
        """
        return self.variant_donor_context is not None

    def get_donor_context_id(self):
        """Return the donor context identifier.

        Returns
        -------
        str
            Donor context identifier
        """
        return self.variant_donor_context.get_context_id()

    def get_donor_context_sample_id(self):
        """Return the sample name/identifier of the donor context.

        Returns
        -------
        str
            Donor context sample name/identifier
        """
        return self.variant_donor_context.get_sample_id()

    def get_donor_context_chrom(self):
        """Return the chromosome of the donor context.

        Returns
        -------
        self.variant_donor_context.get_context_chrom()
            Chromosome of the variant context.
        """
        return self.variant_donor_context.get_context_chrom()

    def get_donor_context_origin(self):
        """Return the origin position of the donor context.

        Returns
        -------
        int
            Variant genomic position that the context is constructed from.
        """
        return self.variant_donor_context.get_context_origin()

    def get_donor_context_start(self):
        """Return the starting position of the donor context.

        Returns
        -------
        int
            Donor context leftmost genomic position
        """
        return self.variant_donor_context.get_context_start()

    def get_donor_context_end(self):
        """Return the ending position of the donor context.

        Returns
        -------
        int
            Donor context rightmost genomic position
        """
        return self.variant_donor_context.get_context_end()

    def get_donor_context_length(self):
        """Return the length of the donor context.

        Returns
        -------
        int
            Donor context length
        """
        return self.variant_donor_context.get_context_length()

    def get_donor_context_reads(self):
        """Return all donor context reads.

        Returns
        -------
        list of pysam.AlignedSegment
            Donor context reads
        """
        return self.variant_donor_context.get_context_bam_reads()

    def get_donor_context_read_ids(self):
        """Return the IDs of all donor context reads.

        Returns
        -------
        list of str
            Donor context read IDs
        """
        return self.variant_donor_context.get_context_read_ids()

    def get_donor_context_read_starts(self):
        """Return the leftmost genomic read positions of all donor context reads.

        Returns
        -------
        list of int
            Donor context leftmost genomic read positions
        """
        return self.variant_donor_context.get_context_read_starts()

    def get_donor_context_read_left_positions(self):
        """Return the leftmost genomic positions of all R1 donor context reads.

        Returns
        -------
        list of int
            Donor context leftmost genomic R1 read positions
        """
        return self.variant_donor_context.get_context_read_left_positions()

    def get_donor_context_read_ends(self):
        """Return the rightmost genomic positions of all donor context reads.

        Returns
        -------
        list of int
            Donor context rightmost genomic read positions
        """
        return self.variant_donor_context.get_context_read_ends()

    def get_donor_context_read_right_positions(self):
        """Return the rightmost genomic positions of all R2 donor context reads.

        Returns
        -------
        list of int
            Donor context rightmost genomic R2 read positions
        """
        return self.variant_donor_context.get_context_read_right_positions()

    def get_donor_context_read_lengths(self):
        """Return the lengths of the donor context reads.

        Returns
        -------
        list of int
            Donor context read lengths
        """
        return self.variant_donor_context.get_context_read_lengths()

    def get_donor_context_unmapped_mate_ids(self):
        """Return the donor context read IDs that have unmapped mates.

        Returns
        -------
        list of str
            Donor context read IDs
        """
        return self.variant_donor_context.get_unmapped_read_mate_ids()

    # ===METHODS TO PRODUCE SOME OUTPUT ABOUT THE VARIANT CONTEXT==============
    def to_string(self):
        """Create the variant context as a String representation.

        The created String representation of the variant context is equal to
        the entry of a variant context file. Each included data attribute is
        separated by a tab.

        Returns
        -------
        str
            Variant context as variant context file entry
        """
        if self.variant_context_areads is None:
            ad_ratio = "N/A"
            list_areads = None
            acon_len = None
            aread_count = 0
        else:
            ad_ratio = float(len(self.variant_context_areads)
                             / len(self.variant_context_dreads))
            areads = list(set(self.get_acceptor_read_ids()))
            areads.sort()
            list_areads = ";".join(areads)
            acon_len = self.variant_acceptor_context.get_context_length()
            aread_count = len(self.variant_context_areads)
        dreads = list(set(self.get_donor_read_ids()))
        dreads.sort()
        list_dreads = ";".join(dreads)
        variants = ";".join([f"{var.chrom}_{var.pos}_{var.ref}_{','.join(var.alts)}"
                             for var in self.variants])
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
                + str(list_dreads) + "\t"
                + variants)

    def to_statistics_string(self):
        """Return a String with basic variant context statistics.

        Returns
        -------
        str
            Basic tab separated variant context statistics
        """
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
