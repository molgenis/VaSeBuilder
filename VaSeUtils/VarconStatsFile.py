import logging
from VarconStats import VarconStats

class VarconStatsFile:
	def __init__(self, statsFileLoc, varconFilter=None):
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
			self.vaseUtilLogger.warning("Could not read varcon stats file " +str(ioe.filename))
		return varconStats
	
	
	# Returns the entire variant context statistics data map
	def getVarconStatsData(self):
		return self.varconStatsData
	
	# Returns an entire VarconStats object containing all statistics for a single variant context id.
	def getVarconStats(self, varconId):
		if(varconId in self.varconStatsData):
			return 
		return None
	
	# Returns the variant context identifier
	def getVariantContextId(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].get_variant_context_id()
		return None
	
	# Returns the average acceptor read length of the variant context
	def getAverageAcceptorReadLength(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorReadLength()
		return None
	
	# Returns the average donor read length of the variant context
	def getAverageDonorReadLength(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorReadLength()
		return None
	
	# Returns the median acceptor read length of the variant context
	def getMedianAcceptorReadLength(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorReadLength()
		return None
	
	# Returns the median donor read length of the variant context
	def getMedianDonorReadLength(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorReadLength()
		return None
	
	# Returns the average acceptor read quality of the variant context
	def getAverageAcceptorReadQuality(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorReadQuality()
		return None
	
	# Returns the average donor read quality of the variant context
	def getAverageDonorReadQuality(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorReadQuality()
		return None
	
	# Returns the median acceptor read quality of the variant context
	def getMedianAcceptorReadQuality(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorReadQuality()
		return None
	
	# Returns the median donor read quality of the variant context
	def getMedianDonorReadQuality(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorReadQuality()
		return None
	
	# Returns the average acceptor mapq values of the variant context
	def getAverageAcceptorMapQ(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageAcceptorMapQ()
		return None
	
	# Returns the average donor mapq values of the variant context
	def getAverageDonorMapQ(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getAverageDonorMapQ()
		return None
	
	# Returns the median acceptor mapq of the variant context
	def getMedianAcceptorMapQ(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianAcceptorMapQ()
		return None
	
	# Returns the median donor mapq of the variant context
	def getMedianDonorMapQ(self, varconId):
		if(varconId in self.varconStatsData):
			return self.varconStatsData[varconId].getMedianDonorMapQ()
		return None
	
	
	# Returns a list of all variant contexts
	def getVariantContextsIds(self):
		return list(self.varconStatsData.keys())
	
	# Returns a list/hashmap of all average acceptor read lengths
	def getAverageAcceptorReadLengths(self, asList=False):
		if(asList):
			return [x.getAverageAcceptorReadLength() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageAcceptorReadLength() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all average donor read lengths
	def getAverageDonorReadLengths(self, asList=False):
		if(asList):
			return [x.getAverageDonorReadLength() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageDonorReadLength() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median ascceptor read lengths
	def getMedianAcceptorReadLengths(self, asList=False):
		if(asList):
			return [x.getMedianAcceptorReadLength() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianAcceptorReadLength() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median donor read lengths
	def getMedianDonorReadLengths(self, asList=False):
		if(asList):
			return [x.getMedianDonorReadLength() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianDonorReadLength() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all average acceptor read qualities
	def getAverageAcceptorReadQualities(self, asList=False):
		if(asList):
			return [x.getAverageAcceptorReadQuality() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageAcceptorReadQuality() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all average donor read qualities
	def getAverageDonorReadQualities(self, asList=False):
		if(asList):
			return [x.getAverageDonorReadQuality() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageDonorReadQuality() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median acceptor read qualities
	def getMedianAcceptorReadQualities(self, asList=False):
		if(asList):
			return [x.getMedianAcceptorReadQuality() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianAcceptorReadQuality() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median donor read qualities
	def getMedianDonorReadQualities(self, asList=False):
		if(asList):
			return [x.getMedianDonorReadQuality() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianDonorReadQuality() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all average acceptor read mapq values
	def getAverageAcceptorMapQs(self, asList=False):
		if(asList):
			return [x.getAverageAcceptorMapQ() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageAcceptorMapQ() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all average donor read mapq values
	def getAverageDonorMapQs(self, asList=False):
		if(asList):
			return [x.getAverageDonorMapQ() for x in list(self.varconStatsData.values())]
		return {k:v.getAverageDonorMapQ() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median acceptor read mapq values
	def getMedianAcceptorMapQs(self, asList=False):
		if(asList):
			return [x.getMedianAcceptorMapQ() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianAcceptorMapQ() for k,v in self.varconStatsData.items()}
	
	# Returns a list/hashmap of all median donor read mapq values
	def getMedianDonorMapQs(self, asList=False):
		if(asList):
			return [x.getMedianDonorMapQ() for x in list(self.varconStatsData.values())]
		return {k:v.getMedianDonorMapQ() for k,v in self.varconStatsData.items()}
