import statistics


class DonorBamRead:
    """The DonorBamRead is used to save data from an aligned read extracted from a BAM or CRAM file.

    DonorBamRead can be used to save the data from aligned reads that were extracted from a BAM or CRAM file.

    Attributes
    ----------
    bam_read_id : str
        Read identifier
    bam_read_flag : int
        The read bitwise flag value
    bam_read_pairnum : str
        Read pair number (1 or 2)
    bam_read_chrom : str
        Chromosome name the read is located on
    bam_read_ref_pos : int
        Leftmost genomic position of the read
    bam_read_length : int
        Read length
    bam_read_end_pos : int
        Read rightmost position (determined by alignment)
    bam_read_cigar : str
        Read CIGAR string
    bam_read_rnext : str
        Chromosome name of the read mate
    bam_read_pnext : int
        Leftmost genomic position of the read mate
    bam_read_tlen : int
        Read TLEN value
    bam_read_seq : str
        Read sequence
    bam_read_qual : str
        Read ASCII quality
    bam_read_map_qual : int
        Read MAPQ value
    """

    def __init__(self, readid, readflag, readpn, readchrom, readstart, readlen, readend, readcigar, readrnext,
                 readpnext, readtlen, readseq, readquals, readmapq):
        """Saves all required data of an aligned read from a BAM or CRAM file.

        Parameters
        ----------
        readid : str
            The read identifier
        readflag : str
            Aligned read bitwise flag
        readpn : str
            The read pair number
        readchrom : str
            The chromosome name the read is located on
        readstart : int
            The leftmost genomic position of the mapped read
        readlen : int
            The total length of the read
        readcigar : str
            Aligned read CIGAR string
        readrnext : str
            Chromsome name the read mate is located on
        readpnext : int
            Leftmost genomic position of the read mate
        readtlen : int
            Aligned read TLEN value
        readseq : str
            The read sequence
        reqdquals : str
            The qualities of the read in ascii
        mapqual : int
            The MAPQ mapping quality
        """
        self.bam_read_id = readid
        self.bam_read_flag = readflag
        self.bam_read_pairnum = readpn
        self.bam_read_chrom = readchrom
        self.bam_read_ref_pos = readstart
        self.bam_read_length = readlen
        self.bam_read_end_pos = readend
        self.bam_read_cigar = readcigar
        self.bam_read_rnext = readrnext
        self.bam_read_pnext = readpnext
        self.bam_read_tlen = readtlen
        self.bam_read_seq = readseq
        self.bam_read_qual = readquals
        self.bam_read_map_qual = readmapq

    # ===METHODS TO GET SAVED DATA FROM THE DONORBAMREAD=======================
    def get_bam_read_id(self):
        """Returns the identifier of the read.

        Returns
        -------
        bam_read_id : str
            The read identifier
        """
        return self.bam_read_id

    def get_bam_read_flag(self):
        """Returns the read alignment flag.

        Returns
        -------
        self.bam_read_flag : str
            Read alignment flag
        """

    def get_bam_read_pair_number(self):
        """Returns the pair number of the read.
        This is 1 for forward, or R1, reads and 2 for reverse, or R2, reads.

        Returns
        -------
        bam_read_pairnum : str
            The read pair number
        """
        return self.bam_read_pairnum

    def get_bam_read_chrom(self):
        """Returns the name of the chromosome where the read has been mapped.

        Returns
        -------
        bam_read_chrom : str
            The chromosome name where the read is mapped
        """
        return self.bam_read_chrom

    def get_bam_read_ref_pos(self):
        """Returns the genomic mapping position (the left most position) of the read.

        Returns
        -------
        bam_read_ref_pos : int
            The read leftmost genomic position
        """
        return self.bam_read_ref_pos

    def get_bam_read_length(self):
        """Returns the length of the read. In VaSeBuilder this length is calculated by means of the CIGAR string by
        using the pysam method infer_read_length()

        Returns
        -------
        bam_read_length : int
            The length of the read
        """
        return self.bam_read_length

    def get_bam_read_ref_end(self):
        """Returns the rightmost genomic position of the read.

        Returns
        -------
        self.bam_read_end_pos : int
            The rightmost position of the read
        """
        return self.bam_read_end_pos

    def get_bam_read_cigar(self):
        """Returns the read alignment CIGAR string

        Returns
        -------
        self.bam_read_cigar : str
            Read alignment CIGAR string
        """
        return self.bam_read_cigar

    def get_bam_read_rnext(self):
        """Returns the chromosome name the read mate is located on.

        Returns
        -------
        self.bam_read_rnext : str
            Chromosome name the read mate is located on
        """
        return self.bam_read_rnext

    def get_bam_read_pnext(self):
        """Returns the leftmost genomic position of the read mate.

        Returns
        -------
        self.bam_read_pnext : ibnt
            Leftmost genomic position of the read mate
        """
        return self.bam_read_pnext

    def get_bam_read_tlen(self):
        """Returns the read alignment TLEN value.

        Returns
        -------
        self.bam_read_tlen : int
            Read alignment TLEN value
        """

    def get_bam_read_sequence(self):
        """Returns the read sequence as a String.

        Returns
        -------
        bam_read_seq : str
            The sequence of the read
        """
        return self.bam_read_seq

    def get_bam_read_qual(self):
        """Returns the read quality.

        Returns
        -------
        bam_read_qual : str
            Read quality
        """
        return self.bam_read_qual

    def get_bam_read_q_scores(self):
        """Converts the String of ASCII quality scores to a list Q-Scores. The Q-Score for each ASCII quality character
        are calculated by obtaining the unicode code point and subtracting 33.

        Returns
        -------
        qscores : list
            Q-Score of the read
        """
        qscores = []
        for qualSymbol in self.bam_read_qual:
            qscores.append(ord(qualSymbol)-33)
        return qscores

    def get_mapping_qual(self):
        """Returns the mapping quality that was assigned to the read.

        Returns
        -------
        bam_read_map_qual : int
            The mapping quality (MAPQ)
        """
        return self.bam_read_map_qual

    def read_has_hardclipped_bases(self):
        """Returns whether the read has hardclipped bases.

        The CIGAR string of the read is checked whether it contains an 'H' (Hardclipped bases).

        Returns
        -------
        bool
            True if read CIGAR contains 'H', False if not
        """
        if self.bam_read_cigar is not None:
            return "H" in self.bam_read_cigar
        return False

    # ===METHOD TO GET STATISTICS DATA FROM THE DONORBAMREAD===================
    def get_average_qscore(self):
        """Calculates and returns the mean Q-Score of the read.

        Returns
        -------
        float
            The mean Q-Score
        """
        qscores = self.get_bam_read_q_scores()
        return statistics.mean(qscores)

    def get_median_qscore(self):
        """Calculates and returns the median Q-Score of the read.

        Returns
        -------
        int
            The median Q-Score of the read
        """
        qscores = self.get_bam_read_q_scores()
        return statistics.median(qscores)

    # ===METHODS TO CHECK WHETHER THE BAM READ IS R1 OR R2=====================
    def is_read1(self):
        """Returns whether the read is the forward/R1 read by checking if the pair number is set to '1'.

        Returns
        -------
        bool
            True if the read has pair number 1, otherwise False
        """
        return self.bam_read_pairnum == "1"

    def is_read2(self):
        """Returns whether the read is the reverse/R2 read by checking if the pair number is set to '2'.

        Returns
        -------
        bool
            True if the read has pair number 2, otherwise False
        """
        return self.bam_read_pairnum == "2"

    # ===METHODS TO RETURN A STRING REPRESENTATION OF THE DONORBAMREAD OBJECT==
    def to_string(self):
        """Returns a String with all the saved data, separated by tabs, of the read.

        Returns
        -------
        str
            String representation of the read.
        """
        return (str(self.bam_read_id) + "\t"
                + str(self.bam_read_pairnum) + "\t"
                + str(self.bam_read_chrom) + "\t"
                + str(self.bam_read_ref_pos) + "\t"
                + str(self.bam_read_length) + "\t"
                + str(self.bam_read_seq) + "\t"
                + str(self.bam_read_qual) + "\t"
                + str(self.bam_read_map_qual))

    def get_as_fastq_seq(self, addpairnum=False):
        """Returns the read as a FastQ entry.

        If the 'addpairnum' argument is set to True the read pair number will be added at the end of the read
        identifier as '/1' or '/2', depending on whether the read is the forward/R1 or reverse/R2 read.

        Returns
        -------
        str
            The read as FastQ entry
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

    def get_as_bam_read(self):
        """Returns the read as a SAM/BAM/CRAM file entry.

        Returns
        -------
        bamcram_entry : str
            Read as SAM/BAM/CRAM entry
        """
        bamcram_entry = f"{self.bam_read_id}\t{self.bam_read_flag}\t{self.bam_read_chrom}\t{self.bam_read_ref_pos}" \
            f"\t{self.bam_read_map_qual}t{self.bam_read_cigar}\t{self.bam_read_rnext}\t{self.bam_read_pnext}" \
            f"\t{self.bam_read_tlen}\t{self.bam_read_seq}\t{self.bam_read_qual}"
        return bamcram_entry
