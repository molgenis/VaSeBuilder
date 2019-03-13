#!/usr/bin/env python
import logging

class GvcfVariant:
	# Saves the variant data.
	def __init__(self, varId, varChr, varPos, varRef, varAlt, varQual, varFilter, varInfo, varFormatSampleData):
		self.variantId = varId	# Variant identifier (can be '.')
		self.variantChrom = varChr
		self.variantPos = varPos
		self.variantRef = varRef
		self.variantAlt = varAlt
		self.variantQual = varQual
		self.variantFilter = varFilter
		self.variantInfo = varInfo	# Hash containing the Field=Value info items
		self.variantFormatSampleData = varFormatSampleData	# Hash containing sampleid=Hash{format_field : sample_data}
		self.calledInfo = {}
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
	
	
	# Returns the identifier of the variant.
	def getVariantId(self):
		return self.variantId
	
	
	# Returns an identifier based on the chromosome and position of the variant (like SNP21_950432). Can be used to give the variant an identifier if the identifier is '.'.
	def getChromPosId(self):
		return "SNP" +str(self.variantChrom)+ "_" +str(self.variantPos)
	
	
	# Returns the chromosome the variant is located on.
	def getVariantChrom(self):
		return self.variantChrom
	
	
	# Returns the chromosomal position the variant is located on.
	def getVariantPos(self):
		return self.variantPos
	
	
	# Returns the reference allele of the variant.
	def getVariantRef(self):
		return self.variantRef
	
	
	# Returns the alternative allele(s) of the variant.
	def getVariantAlt(self):
		return self.variantAlt
	
	
	# Returns the quality score of the variant.
	def getVariantQual(self):
		return self.variantQual
	
	
	# Returns the filter information of the variant.
	def getVariantFilter(self):
		return self.variantFilter
	
	
	# Returns the info fields of the variant.
	def getVariantInfo(self):
		return self.variantInfo
	
	
	# Returns the value of a specified info field.
	def getVariantInfoFieldValue(self, infoField):
		if(self.hasVariantFormatField(infoField)):
			return self.variantInfo[infoField]
		else:
			return None
	
	
	# Returns whether the variant has a certain info field (such as 'AC').
	def hasVariantInfoField(self, infoField):
		if(infoField in self.variantInfo):
			return True
		else:
			return False
	
	
	# Returns the format info of the variant.
	def getVariantFormatSampleData(self):
		return self.variantFormatSampleData
	
	
	# Returns whether a variant.
	def sampleHasFormatField(self, formatfield):
		if(formatfield in self.variantFormatSampleData):
			return True
		else:
			return False
	
	
	# Returns the specified sample information if available, otherwise returns None.
	def getSampleInfo(self, formatfield):
		if(self.sampleHasFormatField(formatfield)):
			return self.variantFormatSampleData[formatfield]
		else:
			return None
	
	
	# Adds the calling information. (Calling information refers to the calling information from the pipeline)
	def setCalledInfo(self, calledInfo):
		self.calledInfo = calledInfo
	
	# Returns whether pipeline called information has been set for this variant.
	def hasCalledInfo(self):
		return bool(self.calledInfo)
	
	# Returns the calling information (information if and whether the variant was called correctly).
	def getCalledInfo(self):
		return self.calledInfo
	
	# Returns the pipeline called genotype of the variant.
	def getCalledGt(self):
		if(self.calledInfo is not None):
			return self.calledInfo['gt']
		return None
	
	# Returns the pipeline called reference allele
	def getCalledRef(self):
		if(self.calledInfo is not None):
			return self.calledInfo['ref']
		return None
	
	# Returns the pipeline called alternative allele
	def getCalledAlt(self):
		if(self.calledInfo is not None):
			return self.calledInfo['alt']
		return None
	
	# Returns the pipeline called position of the variant
	def getCalledPos(self):
		if(self.calledInfo is not None):
			return self.calledInfo['pos']
		return None
	
	# Returns the pipeline calling result (correct/incorrect/not called/indirectly ; cc/ci/nc/ic)
	def getCalledResult(self):
		if(self.calledInfo is not None):
			return self.calledInfo['res']
		return None
	
