"""Pseudo-pysam.AlignedSegment object used as a placeholder when reading varcon files."""


class ReadIdObject:
    """The ReadIdObject is used to save the read identifiers.

    This class is used when a VariantContextFile is read. DonorBamRead can't
    be used as most data is missing.

    Attributes
    ----------
    query_name : str
        Identifier of the read
    pair_number : str
        Pair number of the read
    """

    def __init__(self, readid):
        """Save the read id and set the pairnumber to a default empty value."""
        self.query_name = readid
        self.pair_number = ""

    def get_bam_read_id(self):
        """Return the identifier of the read.

        Returns
        -------
        self.read_id : str
            Identifier of the read
        """
        return self.query_name

    def set_bam_read_id(self, readid):
        """Set the identifier of the read. Overwrite any already existing value.

        Parameters
        ----------
        readid : str
            The read identifier to set
        """
        self.query_name = readid

    def get_bam_read_pair_number(self):
        """Return the read pair number of the read.

        Returns
        -------
        self.pair_number : str
            The read pair number
        """
        return self.pair_number

    def set_bam_read_pair_number(self, pairnum):
        """Set the read pair number of the read. Overwrite any already existing value.

        Parameters
        ----------
        pairnum : str
            The pair number to set
        """
        self.pair_number = pairnum

    def is_read1(self):
        """Return whether the read is the forward/R1 read.

        Returns
        -------
        bool
            True if read pair number is 1., False if not
        """
        return self.pair_number == "1"

    def is_read2(self):
        """Return whether the read is the reverse/R2 read.

        Returns
        -------
        bool
            True if read pair number is 2, False if not
        """
        return self.pair_number == "2"
