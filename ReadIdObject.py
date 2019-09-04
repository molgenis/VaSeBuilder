class ReadIdObject:
    """The ReadIdObject is used to save the read identifiers.

    This class is used when a VariantContextFile is read. DonorBamRead can't be used as most data is missing."""
    def __init__(self, readid):
        """Constructor saves the provided read identifier and sets the pairnumber to a default empty value."""
        self.read_id = readid
        self.pair_number = ""

    # Returns the name/identifier of the read id object
    def get_bam_read_id(self):
        """Returns the identifier of the read."""
        return self.read_id

    # Sets the name/identifier of the read id object
    def set_bam_read_id(self, readid):
        """Sets the identifier of a read."""
        self.read_id = readid

    # Return the pair number of the read id object
    def get_bam_read_pair_number(self):
        """Returns the read pair number of the read."""
        return self.pair_number

    # Sets the pair number of the read id object
    def set_bam_read_pair_number(self, pairnum):
        """Sets the read pair number of the read."""
        self.pair_number = pairnum

    def is_read1(self):
        """Returns whether the read is the forward/R1 read."""
        return self.pair_number == "1"

    def is_read2(self):
        """Returns whether the read is the reverse/R2 read."""
        return self.pair_number == "2"
