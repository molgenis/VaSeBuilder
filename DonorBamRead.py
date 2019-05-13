import statistics

class DonorBamRead:
    # Saves the required BAM read data
    def __init__(self, readId, readPn, readChrom, readStart, readLen, readSeq, readQuals, mapQual):
        self.bamReadId = readId
        self.bamReadPairNum = readPn
        self.bamReadChrom = readChrom
        self.bamReadRefPos = readStart
        self.bamReadLength = readLen
        self.bamReadSeq = readSeq
        self.bamReadQual = readQuals
        self.bamReadMapQual = mapQual



    # ====================METHODS TO GET SAVED DATA FROM THE DONORBAMREAD====================
    # Returns the BAM read identifier.
    def getBamReadId(self):
        return self.bamReadId

    # Returns the BAM read pair number (1 or 2)
    def getBamReadPairNumber(self):
        return self.bamReadPairNum

    # Returns the BAM read chromosome
    def getBamReadChrom(self):
        return self.bamReadChrom

    # Returns the BAM read starting position on the reference sequence
    def getBamReadRefPos(self):
        return self.bamReadRefPos

    # Returns the BAM read length
    def getBamReadLength(self):
        return self.bamReadLength

    # Returns the BAM read ending position on the reference (calculated as starting position + the length of the read)
    def getBamReadRefEnd(self):
        if(self.bamReadLength is not None):
            return (self.bamReadRefPos + self.bamReadLength)
        return -1

    # Returns the BAM read sequence
    def getBamReadSequence(self):
        return self.bamReadSeq

    # Returns the BAM read quality scores.
    def getBamReadQual(self):
        return self.bamReadQual

    # Returns the BAM read quality as an array of Q-Scores.
    def getBamReadQScores(self):
        qscores = []
        for qualSymbol in self.bamReadQual:
            qscores.append(ord(qualSymbol)-33)
        return qscores

    # Returns the maping quality of the BAM read
    def getMappingQual(self):
        return self.bamReadMapQual



    # ====================METHOD TO GET STATISTICS DATA FROM THE DONORBAMREAD====================
    # Returns the average Q-Score
    def getAverageQscore(self):
        qscores = self.getBamReadQScores()
        return statistics.mean(qscores)

    # Returns the median Q-Score
    def getMedianQScore(self):
        qscores = self.getBamReadQScores()
        return statistics.median(qscores)



    # ====================METHODS TO CHECK WHETHER THE BAM READ IS R1 OR R2====================
    # Returns if the BAM read is the first (forward) read
    def isRead1(self):
        return self.bamReadPairNum == '1'

    # Returns if the BAM read is the second (reverse) read
    def isRead2(self):
        return self.bamReadPairNum == '2'



    # ====================METHODS TO RETURN A STRING REPRESENTATION OF THE DONORBAMREAD OBJECT====================
    # Returns a String representation
    def toString(self):
        return str(self.bamReadid) +"\t"+ str(self.bamReadPairNum) +"\t"+ str(self.bamreadChrom) +"\t"+ str(self.bamReadRefPos) +"\t"+ str(self.bamReadLength) +"\t"+ str(self.bamReadSeq) +"\t"+ str(self.bamReadQual) +"\t"+ str(self.bamReadMapQual)

    # Returns the BAM read as a fastq sequence
    def getAsFastQSeq(self, addPairNum=False):
        if(addPairNum):
            return ("@"+ str(self.bamReadId) +"/"+ str(self.bamReadPairNum) +"\n"+ str(self.bamReadSeq) +"\n+\n"+ str(self.bamReadQual) +"\n")
        return ("@"+ str(self.bamReadId) +"\n"+ str(self.bamReadSeq) +"\n+\n"+ str(self.bamReadQual) +"\n")
