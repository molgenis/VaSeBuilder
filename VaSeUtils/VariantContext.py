class VariantContext:
	# Constructor that saves the variant context data.
	def __init__(self, varconType, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconReads, otherVarconReads=None, varconALength=None, varconDLength=None):
		self.variantContextType = varconType
		self.variantContextId = varconId
		self.variantContextSample = varconSample
		self.variantContextChrom = varconChrom
		self.variantContextOrigin = varconOrigin
		self.variantContextStart = varconStart
		self.variantContextEnd = varconEnd
		self.variantContextAcceptorLength = varconALength
		self.variantContextDonorLength = varconDLength
		
		if(varconType=='acceptor'):
			self.variantContextAReadIds = varconReads
			self.variantContextDReadIds = None
			self.variantContextAReadNum = len(varconReads)
			self.variantContextDReadNum = 0
			self.variantContextADRatio = 1.0
		elif(varconType=='donor'):
			self.variantContextAReadIds = None
			self.variantContextDReadIds = varconReads
			self.variantContextAReadNum = 0
			self.variantContextDReadNum = len(varconReads)
			self.variantContextADRatio = 0.0
		else:
			self.variantContextAReadIds = otherVarconReads
			self.variantContextDReadIds = varconReads
			self.variantContextAReadNum = len(otherVarconReads)
			self.variantContextDReadNum = len(varconReads)
			self.variantContextADRatio = float(self.variantContextAReadNum / self.variantContextDReadNum)
		self.variantContextUnmappedMates = None
		self.variantContextStats = None
	
	
	# Returns the type of variant context (acceptor/donor/combined)
	def getVariantContextType(self):
		return self.variantContextType
	
	# Returns the variant context identifier.
	def getVariantContextId(self):
		return self.variantContextId
	
	# Returns the sample identifier the variant context is based on
	def getVariantContextSample(self):
		return self.variantContextSample
	
	# Returns the chromosome name of the variant context
	def getVariantContextChrom(self):
		return self.variantContextChrom
	
	# Returns the position (origin) of the variant the variant context is based on
	def getVariantContextOrigin(self):
		return self.variantContextOrigin
	
	# Returns the starting position of the variant context
	def getVariantContextStart(self):
		return self.variantContextStart
	
	# Returns the ending position of the variant context
	def getVariantContextEnd(self):
		return self.variantContextEnd
	
	# Returns the length of the variant context
	def getVariantContextLength(self):
		return (self.variantContextEnd - self.variantContextStart)
	
	# Returns the length of the acceptor context
	def getAcceptorContextLength(self):
		return self.variantContextAcceptorLength
	
	# Returns the length of the donor context
	def getDonorContextLength(self):
		return self.variantContextDonorLength
	
	# Returns the number of acceptor reads overlapping with the variant context
	def getNumberOfAcceptorReads(self):
		return self.variantContextAReadNum
	
	# Returns the number of donor reads overlapping with the variant context
	def getNumberOfDonorReads(self):
		return self.variantContextDReadNum
	
	# Returns the ratio of acceptor over donor reads for the variant context
	def getAcceptorDonorReadRatio(self):
		return self.variantContextADRatio
	
	# Returns the identifiers of acceptor reads overlapping with the variant context
	def getAcceptorReadIds(self):
		return self.variantContextAReadIds
	
	# Returns the identifiers of donor reads overlapping with the variant context
	def getDonorReadIds(self):
		return self.variantContextDReadIds
	
	# Returns the chormosome and position of the variant the variant context is based on
	def getVariantContextVarPos(self):
		chromPos = self.variantContextId.split('_')
		return [chromPos[0], int(chromPos[1])]
	
	# Returns the varcon statistics object of the current variant context
	def getVariantContextStats(self):
		return self.variantContextStats
	
	
	# Sets the acceptor read id list of the variant context
	def setAcceptorReadIds(self, acceptorReadIds):
		self.variantContextAReadIds = acceptorReadIds
	
	# Sets the donor read id list of the variant context
	def setDonorReadIds(self, donorReadIds):
		self.variantContextDReadIds = donorReadIds
	
	# Sets the length of the acceptor context
	def setAcceptorContextLength(self, aconLength):
		self.variantContextAcceptorLength = aconLength
	
	# Sets the length of the donor context
	def setDonorContextLength(self, dconLength):
		variantContextDonorLength = dconLength
	
	
	# Compares the current VariantContext object to another
	def compare(self, varconObj):
		varconDiff = {}
		# Compare each field of the variant context
		if(self.variantContextId != varconObj.getVariantContextId()):
			varconDiff[1] = [self.variantContextId, varconObj.getVariantContextId()]
		if(self.variantContextSample != varconObj.getVariantContextSample()):
			varconDiff[2] = [self.variantContextSample, varconObj.getVariantContextSample()]
		if(self.variantContextChrom != varconObj.getVariantContextChrom()):
			varconDiff[3] = [self.variantContextChrom, varconObj.getVariantContextChrom()]
		if(self.variantContextOrigin < varconObj.getVariantContextOrigin() or self.variantContextOrigin > varconObj.getVariantContextOrigin()):
			varconDiff[4] = [self.variantContextOrigin, varconObj.getVariantContextOrigin()]
		if(self.variantContextStart < varconObj.getVariantContextStart() or self.variantContextStart > varconObj.getVariantContextStart()):
			varconDiff[5] = [self.variantContextStart, varconObj.getVariantContextStart()]
		if(self.variantContextEnd < varconObj.getVariantContextEnd() or self.variantContextEnd > varconObj.getVariantContextEnd()):
			varconDiff[6] = [self.variantContextEnd, varconObj.getVariantContextEnd()]
		if(self.variantContextAcceptorLength != varconObj.getAcceptorContextLength()):
			varconDiff[7] = [self.variantContextAcceptorLength, varconObj.getAcceptorContextLength()]
		if(self.variantContextDonorLength != varconObj.getDonorContextLength()):
			varconDiff[8] = [self.variantContextDonorLength, varconObj.getDonorContextLength()]
		if(self.variantContextAReadNum != varconObj.getNumberOfAcceptorReads()):
			varconDiff[9] = [self.variantContextAReadNum, varconObj.getNumberOfAcceptorReads()]
		if(self.variantContextDReadNum != varconObj.getNumberOfDonorReads()):
			varconDiff[10] = [self.variantContextDReadNum, varconObj.getNumberOfDonorReads()]
		if(self.variantContextADRatio != varconObj.getAcceptorDonorReadRatio()):
			varconDiff[11] = [self.variantContextADRatio, varconObj.getAcceptorDonorReadRatio()]
		if(self.variantContextAReadIds != varconObj.getAcceptorReadIds()):
			varconDiff[12] = [self.variantContextAReadIds, varconObj.getAcceptorReadIds()]
		if(self.variantContextDReadIds != varconObj.getDonorReadIds()):
			varconDiff[13] = [self.variantContextDReadIds, varconObj.getDonorReadIds()]
		return varconDiff
	
	
	# Returns a String representation 
	def toString(self):
		if(self.variantContextType == 'acceptor'):
			return str(self.variantContextId)+ "\t" +str(self.variantContextSample)+ "\t" +str(self.variantContextChrom)+ "\t" +str(self.variantContextOrigin)+ "\t" +str(self.variantContextStart)+ "\t" +str(self.variantContextEnd)+ "\t" 
			+str(self.variantContextAReadNum)+ "\t" +';'.join(self.variantContextAReadIds)
		if(self.variantContextType == 'donor'):
			return str(self.variantContextId)+ "\t" +str(self.variantContextSample)+ "\t" +str(self.variantContextChrom)+ "\t" +str(self.variantContextOrigin)+ "\t" +str(self.variantContextStart)+ "\t" +str(self.variantContextEnd)+ "\t" 
			+str(self.variantContextDReadNum)+ "\t" +';'.join(self.variantContextDReadIds)
		if(self.variantContextType == 'combi'):
			return str(self.variantContextId)+ "\t" +str(self.variantContextSample)+ "\t" +str(self.variantContextChrom)+ "\t" +str(self.variantContextOrigin)+ "\t" +str(self.variantContextStart)+ "\t" +str(self.variantContextEnd)+ "\t" 
			+str(self.variantContextAcceptorLength)+ "\t" +str(self.variantContextDonorLength)+ "" +str(len(self.variantContextAReadIds))+ "\t" +str(len(self.variantContextDReadIds))+ "\t" 
			+str(float(len(self.variantContextAReadIds)/len(self.variantContextDReadIds)))+ "\t" +';'.join(self.variantContextAReadIds)+ "\t" +';'.join(self.variantContextDReadIds)
		return ""
