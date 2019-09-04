import statistics


class DonorBamRead:
    """The DonorBamRead is used to save data from an aligned read extracted from a BAM or CRAM file.

    DonorBamRead can be used to save  the data from aligned reads that were extracted from a BAM or CRAM file."""
    # Saves the required BAM read data.
    def __init__(self, readid, readpn, readchrom, readstart,
                 readlen, readseq, readquals, mapqual):
        """Saves all required data of an aligned read from a BAM or CRAM file.

        Parameters
        ----------
        readid : str
            The read identifier
        readpn : str
            The read pair number
        readchrom : str
            The chromosome name the read is located on
        readstart : int
            The leftmost genomic position of the mapped read
        readlen : int
            The total length of the read
        readseq : str
            The read sequence
        reqdquals : str
            The qualities of the read in ascii
        mapqual : int
            The MAPQ mapping quality
        """
        self.bam_read_id = readid
        self.bam_read_pairnum = readpn
        self.bam_read_chrom = readchrom
        self.bam_read_ref_pos = readstart
        self.bam_read_length = readlen
        self.bam_read_seq = readseq
        self.bam_read_qual = readquals
        self.bam_read_map_qual = mapqual

    # ===METHODS TO GET SAVED DATA FROM THE DONORBAMREAD=======================
    # Returns the BAM read identifier.
    def get_bam_read_id(self):
        """Returns the identifier of the read.

        Returns
        -------
        str
            The read identifier
        """
        return self.bam_read_id

    # Returns the BAM read pair number (1 or 2).
    def get_bam_read_pair_number(self):
        """Returns the pair number of the read.
        This is 1 for forward, or R1, reads and 2 for reverse, or R2, reads.

        Returns
        -------
        str
            The read pair number
        """
        return self.bam_read_pairnum

    # Returns the BAM read chromosome.
    def get_bam_read_chrom(self):
        """Returns the name of the chromosome where the read has been mapped.

        Returns
        -------
        str
            The chromosome name where the read is mapped
        """
        return self.bam_read_chrom

    # Returns the BAM read starting position on the reference sequence.
    def get_bam_read_ref_pos(self):
        """Returns the genomic mapping position (the left most position) of the read.

        Returns
        -------
        int
            The read leftmost genomic position
        """
        return self.bam_read_ref_pos

    # Returns the BAM read length.
    def get_bam_read_length(self):
        """Returns the length of the read. In VaSeBuilder this length is calculated by means of the CIGAR string by
        using the pysam method infer_read_length()

        Returns
        -------
        int
            The length of the read
        """
        return self.bam_read_length

    # Returns the BAM read ending position on the reference (calculated as starting position + the length of the read).
    def get_bam_read_ref_end(self):
        """Returns the """
        if self.bam_read_length is not None:
            return self.bam_read_ref_pos + self.bam_read_length
        return -1

    # Returns the BAM read sequence.
    def get_bam_read_sequence(self):
        """Returns the read sequence as a String.

        Returns
        -------
        str
            The sequence of the read
        """
        return self.bam_read_seq

    # Returns the BAM read quality scores.
    def get_bam_read_qual(self):
        """Returns the """
        return self.bam_read_qual

    # Returns the BAM read quality as an array of Q-Scores.
    def get_bam_read_q_scores(self):
        """Converts the String of ASCII quality scores to a list Q-Scores. The Q-Score for each ASCII quality character
        are calculated by obtaining the unicode code point and subtracting 33."""
        qscores = []
        for qualSymbol in self.bam_read_qual:
            qscores.append(ord(qualSymbol)-33)
        return qscores

    # Returns the maping quality of the BAM read.
    def get_mapping_qual(self):
        """Returns the mapping quality that was assigned to the read."""
        return self.bam_read_map_qual

    # ===METHOD TO GET STATISTICS DATA FROM THE DONORBAMREAD===================
    # Returns the average Q-Score.
    def get_average_qscore(self):
        """Calculates and returns the mean Q-Score of the read."""
        qscores = self.get_bam_read_q_scores()
        return statistics.mean(qscores)

    # Returns the median Q-Score.
    def get_median_qscore(self):
        """Calculates and returns the median Q-Score of the read."""
        qscores = self.get_bam_read_q_scores()
        return statistics.median(qscores)

    # ===METHODS TO CHECK WHETHER THE BAM READ IS R1 OR R2=====================
    # Returns if the BAM read is the first (forward) read.
    def is_read1(self):
        """Returns whether the read is the forward/R1 read by checking if the pair number is set to '1'."""
        return self.bam_read_pairnum == "1"

    # Returns if the BAM read is the second (reverse) read.
    def is_read2(self):
        """Returns whether the read is the reverse/R2 read by checking if the pair number is set to '2'."""
        return self.bam_read_pairnum == "2"

    # ===METHODS TO RETURN A STRING REPRESENTATION OF THE DONORBAMREAD OBJECT==
    # Returns a String representation.
    def to_string(self):
        """Returns a String with all the saved data, separated by tabs, of the read."""
        return (str(self.bam_read_id) + "\t"
                + str(self.bam_read_pairnum) + "\t"
                + str(self.bam_read_chrom) + "\t"
                + str(self.bam_read_ref_pos) + "\t"
                + str(self.bam_read_length) + "\t"
                + str(self.bam_read_seq) + "\t"
                + str(self.bam_read_qual) + "\t"
                + str(self.bam_read_map_qual))

    # Returns the BAM read as a fastq sequence.
    def get_as_fastq_seq(self, addpairnum=False):
        """Returns the read as a FastQ entry.

        If the 'addpairnum' argument is set to True the read pair number will be added at the end of the read
        identifier as '/1' or '/2', depending on whether the read is the forward/R1 or reverse/R2 read.

        Returns
        -------
        str
            The
        """
        if addpairnum:
            return ("@" + str(self.bam_read_id) + "/" + str(self.bam_read_pairnum) + "\n"
                    + str(self.bam_read_seq) + "\n"
                    + "+\n"
                    + str(self.bam_read_qual) + "\n")

        return ("@" + str(self.bam_read_id) + "\n"
                + str(self.bam_read_seq) + "\n"
                + "+\n"
                + str(self.bam_read_qual) + "\n")
