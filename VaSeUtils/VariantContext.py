class VariantContext:
	# Constructor that saves the variant context data.
	def __init__(self, varconId, varconSample, varconChrom, varconStart, varconEnd):
		self.variantContextId = varconId
		self.variantContextSample = varconSample
		self.variantContextChrom = varconChrom
		self.variantContextStart = varconStart
		self.variantContextEnd = varconEnd
	
	# Returns the variant identifier of the 
	def getVariantContextId(self):
		return self.variantContextId
	
	# Returns the sample identifier the variant context is based on
	def getVariantContextSample(self):
		return self.variantContextSample
	
	# Returns the chromosome name of the variant context
	def getVariantContextChrom(self):
		return self.variantContextChrom
	
	# Returns the starting position of the variant context
	def getVariantContextStart(self):
		return self.variantContextStart
	
	# Returns the ending position of the variant context
	def getVariantContextEnd(self):
		return self.variantContextEnd
	
	# Returns the length of the variant context
	def getVariantContextLength(self):
		return (self.variantContextEnd - self.variantContextStart)
	
	# Returns the chormosome and position of the variant the variant context is based on
	def getVariantContextVarPos(self):
		if(self.variantContextId.startswith('SNP')):
			chromPos = self.variantContextId.split('_')
			return [chromPos[0][3:], int(chromPos[1])]
		return ['-1', -1]
