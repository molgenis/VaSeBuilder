import statistics


class DonorBamRead:
    # Saves the required BAM read data.
    def __init__(self, readid, readpn, readchrom, readstart,
                 readlen, readseq, readquals, mapqual):
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
        return self.bam_read_id

    # Returns the BAM read pair number (1 or 2).
    def get_bam_read_pair_number(self):
        return self.bam_read_pairnum

    # Returns the BAM read chromosome.
    def get_bam_read_chrom(self):
        return self.bam_read_chrom

    # Returns the BAM read starting position on the reference sequence.
    def get_bam_read_ref_pos(self):
        return self.bam_read_ref_pos

    # Returns the BAM read length.
    def get_bam_read_length(self):
        return self.bam_read_length

    # Returns the BAM read ending position on the reference (calculated
    # as starting position + the length of the read).
    def get_bam_read_ref_end(self):
        if self.bam_read_length is not None:
            return self.bam_read_ref_pos + self.bam_read_length
        return -1

    # Returns the BAM read sequence.
    def get_bam_read_sequence(self):
        return self.bam_read_seq

    # Returns the BAM read quality scores.
    def get_bam_read_qual(self):
        return self.bam_read_qual

    # Returns the BAM read quality as an array of Q-Scores.
    def get_bam_read_q_scores(self):
        qscores = []
        for qualSymbol in self.bam_read_qual:
            qscores.append(ord(qualSymbol)-33)
        return qscores

    # Returns the maping quality of the BAM read.
    def get_mapping_qual(self):
        return self.bam_read_map_qual

    # ===METHOD TO GET STATISTICS DATA FROM THE DONORBAMREAD===================
    # Returns the average Q-Score.
    def get_average_qscore(self):
        qscores = self.get_bam_read_q_scores()
        return statistics.mean(qscores)

    # Returns the median Q-Score.
    def get_median_qscore(self):
        qscores = self.get_bam_read_q_scores()
        return statistics.median(qscores)

    # ===METHODS TO CHECK WHETHER THE BAM READ IS R1 OR R2=====================
    # Returns if the BAM read is the first (forward) read.
    def is_read1(self):
        return self.bam_read_pairnum == "1"

    # Returns if the BAM read is the second (reverse) read.
    def is_read2(self):
        return self.bam_read_pairnum == "2"

    # ===METHODS TO RETURN A STRING REPRESENTATION OF THE DONORBAMREAD OBJECT==
    # Returns a String representation.
    def to_string(self):

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
        if addpairnum:
            return ("@" + str(self.bam_read_id) + "/" + str(self.bam_read_pairnum) + "\n"
                    + str(self.bam_read_seq) + "\n"
                    + "+\n"
                    + str(self.bam_read_qual) + "\n")

        return ("@" + str(self.bam_read_id) + "\n"
                + str(self.bam_read_seq) + "\n"
                + "+\n"
                + str(self.bam_read_qual) + "\n")
