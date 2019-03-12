#!/use/bin/env python
import logging
from GvcfVariant import GvcfVariant

class GvcfFile:
	# Constructor reading and saving 
	def __init__(self, gvcfFileLoc):
		self.sampleIds = []
		self.gvcfData = {}
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
		self.readGvcfFile(gvcfFileLoc)
	
	
	# Reads the GVCF file produced by the NGS_DNA pipeline and saves the results.
	def readGvcfFile(self, gvcfFileLoc):
		try:
			with open(gvcfFileLoc, 'r') as vcfFile:
				for vcfLine in vcfFile:
					if(not vcfLine.startswith("##")):
						vcfLine = vcfLine.strip()
						vcfLineData = vcfLine.split("\t")
						if(vcfLine.startswith("#CHROM")):
							self.sampleIds = vcfLineData[9:]	# Extract the variant info
						else:
							# Save the data in a GvcfVariant object
							varInfoData = self.getVariantInfo(vcfLineData[7])
							varFormatSampleData = self.getVariantFormatSampleData(vcfLineData[8], vcfLineData[9])
							variantId = self.getVariantId(vcfLineData[2], vcfLineData[0], vcfLineData[1])
							self.gvcfData[variantId] = GvcfVariant(vcfLineData[2], vcfLineData[0], vcfLineData[1], vcfLineData[3], vcfLineData[4], vcfLineData[5], vcfLineData[6], varInfoData, varFormatSampleData)
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not read vcf/gvcf file " +gvcfFileLoc)
			exit()
	
	
	# Returns an identifier for the SNP variant.
	def getVariantId(self, varId, varChrom, varPos):
		if(varId=='.'):
			return "SNP" +str(varChrom)+ "_" +str(varPos)
		return varId
	
	
	# Returns the info within the INFO column of the variant
	def getVariantInfo(self, variantInfoData):
		infoData = {}
		varInfoPairs = variantInfoData.split(';')
		for varInfoPair in varInfoPairs:
			infoField, infoValue = varInfoPair.split('=')
			infoData[infoField] = float(infoValue)
		return infoData
	
	
	# Returns the data within the format column and sample column(s).
	def getVariantFormatSampleData(self, varFormatFields, varSampleFields):
		formatSampleData = {}
		formatData = varFormatFields.split(':')	# Split the format and sample data on ':'
		sampleDataElements = varSampleFields.split(':')
			
		# Add the format, sample data fields for the current sample.
		for fieldIndex in range(0, len(formatData)):
			formatSampleData[formatData[fieldIndex]] = sampleDataElements[fieldIndex]
		return formatSampleData
	
	
	# Returns the number of read variants
	def getNumberOfVariants(self):
		return len(self.gvcfData)
	
	
	# Returns the list of sample ids within the gVCF.
	def getSampleIds(self):
		return self.sampleIds
	
	
	# Returns the variant data saved in the file
	def getVariantData(self):
		return self.gvcfData
	
	
	
	# Checks whether the variants of one file are present in the other and have the same genotype
	def checkVariantsInOther(self, resultsGvcf, variantContextFile):
		callingInfo = {'cc' : [], 'ci' : [], 'nc' : [], 'ic' : []}	# cc=called correctly, ci=called incorrectly, nc=not called, ic=indirectly called
		resultsGvcfData = resultsGvcf.getVariantData()
		
		for vcfVariantId in self.gvcfData:
			# Check if the variant is in the pipeline result GVCF data.
			resultsVariantid = self.getVariantIdentifierInOther(self.gvcfData[vcfVariantId].getVariantId(), self.gvcfData[vcfVariantId].getChromPosId(), resultsGvcfData)
			callResult = ''
			
			# Determine what type of call it is (cc/ci/nc/ic)
			if(variantid is not None):
				# Check if the variant is called correctly.
				if(variantContextFile.hasVariantContext(vcfVariantId)):
					if(self.variantCalledCorrectly(resultsVariantid, self.gvcfData[vcfVariantId], resultsGvcfData)):
						callingInfo['cc'].append(vcfVariantId)
						callResult = 'cc'
					else:
						callingInfo['ci'].append(vcfVariantId)
						callResult = 'ci'
				else:
					callingInfo['ic'].append(vcfVariantId)
					callResult = 'ic'
				
				# Set the calling information for this particular variant and add it to the donor variant.
				calledInformation = {}
				calledInformation['gt'] = resultsGvcfData[variantid].getSampleInfo('GT')
				calledInformation['ref'] = resultsGvcfData[variantid].getVariantRef()
				calledInformation['alt'] = resultsGvcfData[variantid].getVariantAlt()
				calledInformation['pos'] = resultsGvcfData[variantid].getVariantPos()
				calledInformation['res'] = callResult
				self.gvcfData[vcfVariantId].setCalledInfo(calledInformation)
			else:
				callingInfo['nc'].append(vcfVariantId)
		return callingInfo
	
	
	# Returns which variant identifier is within the other GvcfFile object data. Returns None if not present.
	def getVariantIdentifierInOther(self, variantId, variantChromPosId, gvcfVariantData):
		if(variantId in gvcfVariantData):
			return variantId
		elif(variantChromPosId in gvcfVariantData):
			return variantChromPosId
		else:
			return None
	
	
	# Checks (REF/ALT/GT) whether a specified variant has been called correctly in the other data.
	def variantCalledCorrectly(self, variantId, vcfVariant, gvcfVariantData):
		refCheck = vcfVariant.getVariantRef() == gvcfVariantData[variantId].getVariantRef()
		altCheck = vcfVariant.getVariantAlt() == gvcfVariantData[variantId].getVariantAlt()
		gtCheck = vcfVariant.getSampleInfo('GT') == gvcfVariantData[variantId].getSampleInfo('GT')
		
		#Perform the three checks
		if(refCheck and altCheck and gtCheck):
			return True
		return False
	
	
	# Returns whether a variant is located in the VCF file.
	def hasVariant(self, variantId, variantChrom, variantPos):
		if(variantId in self.gvcfData):
			return True
		elif("SNP"+str(variantChrom)+"_"+str(variantPos) in self.gvcfData):
			return True
		else:
			return False
	
	
	# Returns the position of a specified variant
	def getVariantPos(self, varId):
		if(varId in self.gvcfData):
			return self.gvcfData[varId].getVariantPos()
		return None
	
	# Returns the position of a variant if present in the file.
	def getVariantPos2(self, varId, varChr, varPos):
		variantid = self.getVariantId(varId, varChr, varPos)
		if(variantid in self.gvcfData):
			return self.gvcfData[variantid].getVariantPos()
		return None
	
	# Returns the reference allele of a specified variant
	def getVariantRef(self, varId):
		if(varId in self.gvcfData):
			return self.gvcfData[varId].getVariantRef()
		return None
	
	# Returns the reference allele of a variant if present in the file.
	def getVariantRef2(self, varId, varChr, varPos):
		variantid = self.getVariantId(varId, varChr, varPos)
		if(variantid in self.gvcfData):
			return self.gvcfData[variantid].getVariantRef()
		return None
	
	# Returns the alternative allele of a specified variant
	def getVariantAlt(self, varId):
		if(varId in self.gvcfData):
			return self.gvcfData[varId].getVariantAlt()
		return None
	
	# Returs the alternative allele of the variant if present in the file.
	def getVariantAlt2(self, varId, varChr, varPos):
		variantid = self.getVariantId(varId, varChr, varPos)
			if(variantid in self.gvcfData):
				return self.gvcfData[variantid].getVariantAlt()
			return None
	
	# Returns the quality of a specified variant
	def getVariantQual(self, varId):
		if(varId in self.gvcfData):
			return self.gvcfData[varId].getVariantQual()
		return None
	
	# Returns the variant quality if present in the file.
	def getVariantQual2(self, varId, varChr, varPos):
		variantid = self.getVariantId(varId, varChr, varPos)
		if(variantid in self.gvcfData):
			return self.gvcfData[variantid].getVariantQual()
		return None
	
	# Returns the variant filter info of a specified variant.
	def getVariantFilter(self, varId):
		if(varId in self.gvcfData):
			return self.gvcfData[varId]
		return None
	
	# Returns the variant filter information if present in the file.
	def getVariantFilter2(self, varId, varChr, varPos):
		variantid = self.getVariantId(varId, varChr, varPos)
		if(variantid in self.gvcfData):
			return self.gvcfData[variantid].getVariantFilter()
		return None
