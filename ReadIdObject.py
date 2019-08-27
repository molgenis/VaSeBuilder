class ReadIdObject:
    def __init__(self, readid):
        self.read_id = readid
        self.pair_number = ""

    # Returns the name/identifier of the read id object
    def get_bam_read_id(self):
        return self.read_id

    # Sets the name/identifier of the read id object
    def set_bam_read_id(self, readid):
        self.read_id = readid

    # Return the pair number of the read id object
    def get_bam_read_pair_number(self):
        return self.pair_number

    # Sets the pair number of the read id object
    def set_bam_read_pair_number(self, pairnum):
        self.pair_number = pairnum
