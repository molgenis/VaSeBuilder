class VariantContext:
	# Constructor that saves the variant context data.
	def __init__(self, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconALength, varconDLength, varconAnum, varconDnum, varconAdRatio, varconAReads, varconDReads):
		self.variantContextId = varconId
		self.variantContextSample = varconSample
		self.variantContextChrom = varconChrom
		self.variantContextOrigin = varconOrigin
		self.variantContextStart = varconStart
		self.variantContextEnd = varconEnd
		self.variantContextAcceptorLength = varconALength
		self.variantContextDonorLength = varconDLength
		self.variantContextAReadNum = varconAnum
		self.variantContextDReadNum = varconDnum
		self.variantContextADRatio = varconAdRatio
		self.variantContextAReadIds = varconAReads
		self.variantContextDReadIds = varconDReads
	
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
	
	# Returns the identifiers of donor reads overlapping with the variantg context
	def getDonorReadIds(self):
		return self.variantContextDReadIds
	
	# Returns the chormosome and position of the variant the variant context is based on
	def getVariantContextVarPos(self):
		chromPos = self.variantContextId.split('_')
		return [chromPos[0][3:], int(chromPos[1])]
	
	# Compares the current VariantContext object to another
	def compare(self, varconObj):
		varconDiff = []
		# Compare each field of the variant context
		if(self.variantContextId != varconObj.getVariantContextId()):
			varconDiff.append(1)
		if(self.variantContextSample != varconObj.getVariantContextSample()):
			varconDiff.append(2)
		if(self.variantContextChrom != varconObj.getVariantContextChrom()):
			varconDiff.append(3)
		if(self.variantContextOrigin < varconObj.getVariantContextOrigin() or self.variantContextOrigin > varconObj.getVariantContextOrigin()):
			varconDiff.append(4)
		if(self.variantContextStart < varconObj.getVariantContextStart() or self.variantContextStart > varconObj.getVariantContextStart()):
			varconDiff.append(5)
		if(self.variantContextEnd < varconObj.getVariantContextEnd() or self.variantContextEnd > varconObj.getVariantContextEnd()):
			varconDiff.append(6)
		if(self.variantContextAcceptorLength != varconObj.getAcceptorContextLength()):
			varconDiff.append(7)
		if(self.variantContextDonorLength ! = varconObj.getDonorContextLength()):
			varconDiff.append(8)
		if(self.variantContextAReadNum != varconObj.getNumberOfAcceptorReads()):
			varconDiff.append(9)
		if(self.variantContextDReadNum != varconObj.getNumberOfDonorReads()):
			varconDiff.append(10)
		if(self.variantContextADRatio != varconObj.getAcceptorDonorReadRatio()):
			varconDiff.append(11)
		if(self.variantContextAReadIds != varconObj.getAcceptorReadIds()):
			varconDiff.append(12)
		if(self.variantContextDReadIds != varconObj.getDonorReadIds()):
			varconDiff.append(13)
		return varconDiff
