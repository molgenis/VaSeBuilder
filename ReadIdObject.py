class ReadIdObject:
    """The ReadIdObject is used to save the read identifiers.

    This class is used when a VariantContextFile is read. DonorBamRead can't be used as most data is missing.

    Attributes
    ----------
    read_id : str
        Identifier of the read
    pair_number : str
        Pair number of the read
    """
    def __init__(self, readid):
        """Constructor saves the provided read identifier and sets the pairnumber to a default empty value."""
        self.read_id = readid
        self.pair_number = ""

    def get_bam_read_id(self):
        """Returns the identifier of the read.

        Returns
        -------
        self.read_id : str
            Identifier of the read
        """
        return self.read_id

    def set_bam_read_id(self, readid):
        """Sets the identifier of the read. Overwrites any already existing value.

        Parameters
        ----------
        readid : str
            The read identifier to set
        """
        self.read_id = readid

    def get_bam_read_pair_number(self):
        """Returns the read pair number of the read.

        Returns
        -------
        self.pair_number : str
            The read pair number
        """
        return self.pair_number

    def set_bam_read_pair_number(self, pairnum):
        """Sets the read pair number of the read. Overwrites any already existing value.

        Parameters
        ----------
        pairnum : str
            The pair number to set
        """
        self.pair_number = pairnum

    def is_read1(self):
        """Returns whether the read is the forward/R1 read.

        Returns
        -------
        bool
            True if read pair number is 1., False if not
        """
        return self.pair_number == "1"

    def is_read2(self):
        """Returns whether the read is the reverse/R2 read.

        Returns
        -------
        bool
            True if read pair number is 2, Fale if not
        """
        return self.pair_number == "2"
