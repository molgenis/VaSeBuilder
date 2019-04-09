class VariantContext:
	# Constructor that saves the variant context data.
	def __init__(self, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd):
		self.variantContextId = varconId
		self.variantContextSample = varconSample
		self.variantContextChrom = varconChrom
		self.variantContextOrigin = varconOrigin
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
	
	# Returns the chormosome and position of the variant the variant context is based on
	def getVariantContextVarPos(self):
		if(self.variantContextId.startswith('SNP')):
			chromPos = self.variantContextId.split('_')
			return [chromPos[0][3:], int(chromPos[1])]
		return ['-1', -1]
	
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
		return varconDiff
