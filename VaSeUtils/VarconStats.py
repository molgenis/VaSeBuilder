class VarconStats:
	def __init__(self, varconId, avgArLength, avgDrLength, medArLength, medDrLength, avgArQual, avgDrQual, medArQual, medDrQual, avgAMapQ, avgDMapQ, medAMapQ, medDMapQ):
		self.variantContextId = varconId
		self.varconAvgAReadLength = avgArLength
		self.varconAvgDReadLength = avgDrLength
		self.varconMedianAReadLength = medArLength
		self.varconMedianDReadLength = medDrLength
		self.varconAvgAReadQual = avgArQual
		self.varconAvgDReadQual = avgDrQual
		self.varconMedianAReadQual = medArQual
		self.varconMedianDReadQual = medDrQual
		self.varconAvgAMapQ = avgAMapQ
		self.varconAvgDMapQ = avgDMapQ
		self.varconMedianAMapQ = medAMapQ
		self.varconMedianDMapQ = avgDMapQ
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	# Returns the variant context identifier
	def getVariantContextId(self):
		return self.variantContextId
	
	# Returns the average acceptor read length of the variant context
	def getAverageAcceptorReadLength(self):
		return self.varconAvgAReadLength
	
	# Returns the average donor read length of the variant context
	def getAverageDonorReadLength(self):
		return self.varconAvgDReadLength
	
	# Returns the median acceptor read length of the variant context
	def getMedianAcceptorReadLength(self):
		return self.varconMedianAReadLength
	
	# Returns the median donor read length of the variant context
	def getMedianDonorReadLength(self):
		return self.varconMedianDReadLength
	
	# Returns the average acceptor read quality of the variant context
	def getAverageAcceptorReadQuality(self):
		return self.varconAvgAReadQual
	
	# Returns the average donor read quality of the variant context
	def getAverageDonorReadQuality(self):
		return self.varconAvgDReadQual
	
	# Returns the median acceptor read quality of the variant context
	def getMedianAcceptorReadQuality(self):
		return self.varconMedianAReadQUal
	
	# Returns the median donor read quality of the variant context
	def getMedianDonorReadQuality(self):
		return self.varconMedianDReadQual
	
	# Returns the average acceptor mapq values of the variant context
	def getAverageAcceptorMapQ(self):
		return self.varconAvgAMapQ
	
	# Returns the average donor mapq values of the variant context
	def getAverageDonorMapQ(self):
		return self.varconAvgDMapQ
	
	# Returns the median acceptor mapq of the variant context
	def getMedianAcceptorMapQ(self):
		return self.varconMedianAMapQ
	
	# Returns the median donor mapq of the variant context
	def getMedianDonorMapQ(self):
		return self.varconMedianDMapQ
	
	# Returns the variant context statistics as a file entry
	def tostring(self):
		return str(self.variantContextId)+ "\t" +str(varconAvgAReadLength)+ "\t" +str(varconAvgDReadLength)+ "" +str()+ "" +str()+ "" +str()+ "" +str()+ "" +str()+ "" +str()+ "" +str()+ ""
