import logging

class VarconStatsFile:
	def __init__(self, statsFileLoc):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.varconStatsData = self.readVarconStatsFile(statsFileLoc)
	
	
	# Reads the variant context statistics file.
	def readVarconStatsFile(self, fileLoc):
		varconStats = {}
		try:
			with open(fileLoc, 'r') as varconStatsFile:
				next(varconStatsFile)	# Skips the header line
				for fileLine in varconStatsFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					varconStats[fileLineData[0]] = VarconStats(fileLineData[0], fileLineData[1], fileLineData[2], fileLineData[3], fileLineData[4], fileLineData[5], fileLineData[6], fileLineData[7], fileLineData[8], fileLineData[9], fileLineData[10], fileLineData[11], fileLineData[12])
		except IOError as ioe:
			self.vaseUtilLogger.warning("Could not read ")
		return varconStats
	
	
	# Returns the entire variant context statistics data map
	def getVarconStatsData(self):
		return self.varconStatsData
	
	def getVarconStats(self, varconId):
		if(varconId in self.varconStatsData):
			return 
		return None
	
	# Returns the variant context identifier
	def getVariantContextId(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getVariantContextId()
		return None
	
	# Returns the average acceptor read length of the variant context
	def getAverageAcceptorReadLength(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorReadLength()
		return None
	
	# Returns the average donor read length of the variant context
	def getAverageDonorReadLength(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorReadLength()
		return None
	
	# Returns the median acceptor read length of the variant context
	def getMedianAcceptorReadLength(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorReadLength()
		return None
	
	# Returns the median donor read length of the variant context
	def getMedianDonorReadLength(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorReadLength()
		return None
	
	# Returns the average acceptor read quality of the variant context
	def getAverageAcceptorReadQuality(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorReadQuality()
		return None
	
	# Returns the average donor read quality of the variant context
	def getAverageDonorReadQuality(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorReadQuality()
		return None
	
	# Returns the median acceptor read quality of the variant context
	def getMedianAcceptorReadQuality(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorReadQuality()
		return None
	
	# Returns the median donor read quality of the variant context
	def getMedianDonorReadQuality(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorReadQuality()
		return None
	
	# Returns the average acceptor mapq values of the variant context
	def getAverageAcceptorMapQ(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorMapQ()
		return None
	
	# Returns the average donor mapq values of the variant context
	def getAverageDonorMapQ(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorMapQ()
		return None
	
	# Returns the median acceptor mapq of the variant context
	def getMedianAcceptorMapQ(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorMapQ()
		return None
	
	# Returns the median donor mapq of the variant context
	def getMedianDonorMapQ(self):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorMapQ()
		return None
