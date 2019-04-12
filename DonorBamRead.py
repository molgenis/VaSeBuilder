class DonorBamRead:
	def __init__(self, readId, readPn, readChrom, readStart, readLen, readSeq, readQuals):
		self.bamReadId = readId
		self.bamReadPairNum = readPn
		self.bamReadChrom = readChrom
		self.bamReadRefPos = readStart
		self.bamReadLength = readLen
		self.bamReadSeq = readSeq
		self.bamReadQual = readQuals
	
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
	
	# Returns if the BAM read is the first (forward) read
	def isRead1(self):
		return self.bamReadPairNum == '1'
	
	# Returns if the BAM read is the second (reverse) read
	def isRead2(self):
		return self.bamReadPairNum == '2'
	
	# Returns the BAM read as a fastq sequence
	def getAsFastQSeq(self, addPairNum=False):
		if(addPairNum):
			return ("@"+ self.bamReadId +"/"+ self.bamReadPairNum +"\n"+ self.bamReadSeq +"\n+\n"+ self.bamReadQual +"\n")
		return ("@"+ self.bamReadId +"\n"+ self.bamReadSeq +"\n+\n"+ self.bamReadQual +"\n")
