import statistics


class DonorBamRead:
    # Saves the required BAM read data.
    def __init__(self, readId, readPn, readChrom, readStart,
                 readLen, readSeq, readQuals, mapQual):
        self.bamReadId = readId
        self.bamReadPairNum = readPn
        self.bamReadChrom = readChrom
        self.bamReadRefPos = readStart
        self.bamReadLength = readLen
        self.bamReadSeq = readSeq
        self.bamReadQual = readQuals
        self.bamReadMapQual = mapQual

    # ===METHODS TO GET SAVED DATA FROM THE DONORBAMREAD=======================
    # Returns the BAM read identifier.
    def get_bam_read_id(self):
        return self.bamReadId

    # Returns the BAM read pair number (1 or 2).
    def get_bam_read_pair_number(self):
        return self.bamReadPairNum

    # Returns the BAM read chromosome.
    def get_bam_read_chrom(self):
        return self.bamReadChrom

    # Returns the BAM read starting position on the reference sequence.
    def get_bam_read_ref_pos(self):
        return self.bamReadRefPos

    # Returns the BAM read length.
    def get_bam_read_length(self):
        return self.bamReadLength

    # Returns the BAM read ending position on the reference (calculated
    # as starting position + the length of the read).
    def get_bam_read_ref_end(self):
        if (self.bamReadLength is not None):
            return (self.bamReadRefPos + self.bamReadLength)
        return -1

    # Returns the BAM read sequence.
    def get_bam_read_sequence(self):
        return self.bamReadSeq

    # Returns the BAM read quality scores.
    def get_bam_read_qual(self):
        return self.bamReadQual

    # Returns the BAM read quality as an array of Q-Scores.
    def get_bam_read_q_scores(self):
        qscores = []
        for qualSymbol in self.bamReadQual:
            qscores.append(ord(qualSymbol)-33)
        return qscores

    # Returns the maping quality of the BAM read.
    def get_mapping_qual(self):
        return self.bamReadMapQual

    # ===METHOD TO GET STATISTICS DATA FROM THE DONORBAMREAD===================
    # Returns the average Q-Score.
    def get_average_qscore(self):
        qscores = self.get_bam_read_q_scores()
        return statistics.mean(qscores)

    # Returns the median Q-Score.
    def get_median_q_score(self):
        qscores = self.get_bam_read_q_scores()
        return statistics.median(qscores)

    # ===METHODS TO CHECK WHETHER THE BAM READ IS R1 OR R2=====================
    # Returns if the BAM read is the first (forward) read.
    def is_read1(self):
        return self.bamReadPairNum == '1'

    # Returns if the BAM read is the second (reverse) read.
    def is_read2(self):
        return self.bamReadPairNum == '2'

    # ====METHODS TO RETURN A STRING REPRESENTATION OF THE DONORBAMREAD OBJECT=
    # Returns a String representation.
    def to_string(self):

        return (str(self.bamReadId) + "\t"
                + str(self.bamReadPairNum) + "\t"
                + str(self.bamReadChrom) + "\t"
                + str(self.bamReadRefPos) + "\t"
                + str(self.bamReadLength) + "\t"
                + str(self.bamReadSeq) + "\t"
                + str(self.bamReadQual) + "\t"
                + str(self.bamReadMapQual))

    # Returns the BAM read as a fastq sequence.
    def get_as_fastq_seq(self, addPairNum=False):
        if (addPairNum):
            return ("@" + str(self.bamReadId) + "/" + str(self.bamReadPairNum)
                    + "\n"
                    + str(self.bamReadSeq)
                    + "\n+\n"
                    + str(self.bamReadQual)
                    + "\n")
        return ("@" + str(self.bamReadId)
                + "\n"
                + str(self.bamReadSeq)
                + "\n+\n"
                + str(self.bamReadQual)
                + "\n")
