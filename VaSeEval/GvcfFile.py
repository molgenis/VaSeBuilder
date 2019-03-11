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
		
		for vcfVariant in self.gvcfData:
			# Check if the variant is in the gvcf data.
			variantid = self.getVariantIdentifierInOther(vcfVariant.getVariantId(), vcfVariant.getChromPosId(), resultsGvcfData)
			if(variantid is not None):
				# Check if the variant is called correctly.
				if(variantContextFile.hasVariantContext(variantid)):
					if(self.variantCalledCorrectly(variantid, vcfVariant, resultsGvcfData)):
						callingInfo['cc'].append(variantid)
					else:
						callingInfo['ci'].append(variantid)
				else:
					callingInfo['ic'].append(variantid)
				
				# Set the calling information for this particular variant.
				calledInformation = {}
				calledInformation['gt'] = resultsGvcfData[variantid].getSampleInfo('GT')
				calledInformation['ref'] = resultsGvcfData[variantid].getVariantRef()
				calledInformation['alt'] = resultsGvcfData[variantid].getVariantAlt()
				calledInformation['pos'] = resultsGvcfData[variantid].getVariantPos()
				vcfVariant.setCalledInfo(calledInformation)
			else:
				# Do something with variants not found by the pipeline
				if(vcfVariant.getVariantId() == '.'):
					callingInfo['nc'].append(vcfVariant.getChromPosId())
				else:
					callingInfo['nc'].append(vcfVariant.getVariantId())
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
