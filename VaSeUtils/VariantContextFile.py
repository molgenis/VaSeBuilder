import logging
from VariantContext import VariantContext

# Class that can be used to read and construct a variant context file
class VariantContextFile:
	def __init__(self, varconType, fileLoc=None, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.variantContextFileType = varconType
		self.variantContextFileLocation = fileLoc
		self.variantContextsBySample = {}
		self.variantContextsById = {}
		self.variantContexts = []
		self.variantContextStatistics = None
		self.varconFields = {1: 'variant context id',
			2 : 'sample id',
			3 : 'chromosome',
			4 : 'origin',
			5 : 'start pos',
			6 : 'end pos',
			7 : 'acceptor context',
			8 : 'donor context',
			9 : 'number of acceptor reads',
			10 : 'number of donor reads',
			11 : 'acceptor/donor ratio',
			12 : 'acceptor read ids',
			13 : 'donor read ids'
			}
		if(fileLoc is not None):
			self.readVariantContextFile(self, fileLoc, sampleFilter, varconFilter, chromFilter, posFilter)	# Read the provided variant context file with set optional filters
	
	
	# Reads a provided variant context file and saves data according to set filters
	def readVariantContextFile(self, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None):
		try:
			with open(fileLoc, 'r') as vcFile:
				next(vcFile)	# Skip the header line
				for fileLine in vcFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.passesFilter(fileLineData[0], varconFilter)
					chromPass = self.passesFilter(fileLineData[2], chromFilter)
					
					if(samplePass and varconPass and chromPass):
						varconObj = None
						if(self.variantContextFileType=='acceptor' or self.variantContextFileType=='donor'):
							varconObj = VariantContext(self.variantContextFileType, fileLineData[0], fileLineData[1], fileLineData[2], int(fileLineData[3]), int(fileLineData[4]), int(fileLineData[5]), fileLineData[7].split(';'))
						else:
							varconObj = varconObj = VariantContext(self.variantContextFileType, fileLineData[0], fileLineData[1], fileLineData[2], int(fileLineData[3]), int(fileLineData[4]), int(fileLineData[5]), int(fileLineData[6]), int(fileLineData[7]), int(fileLineData[8]), int(fileLineData[9]), float(fileLineData[10]), fileLineData[11].split(';'), fileLineData[12].split(';'))
						if(fileLineData[1] not in self.variantContextsBySample):
							self.variantContextsBySample[fileLineData[1]] = []
						self.variantContextsBySample[fileLineData[1]].append(varconObj)
						self.variantContextsById[fileLineData[0]] = varconObj
						self.variantContexts.append(varconObj)
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read varcon file " +str(ioe.filename))
			exit()
	
	# Compares a set or all variant contexts of the current Variant Context file to another and returns the differences
	def compare(self, varconFile, varconIds=[]):
		vcFileDiffs = {}
		otherVarcon = varconFile.getVariantContextsById()
		
		# Start comparing the two files completely or specific set of variant contexts in both files
		if(len(varconIds)>0):
			for varconId in varconIds:
				vcdiffs = self.variantContextsById[varconId].compare(otherVarcon[varconId])
				if(len(vcdiffs)>0):
					vcFileDiffs[varconId] = vcdiffs
		else:
			for varconId, varconObj in self.variantContextsById.items():
				vcdiffs = varconObj.compare(otherVarcon[varconId])
				if(len(vcdiffs)>0):
					vcFileDiffs[varconId] = vcdiffs
		return vcfFileDiffs
	
	
	# Returns the variant position based on the chr_pos identifier
	def getVariantContextVarPos(self):
		chromPos = self.variantContextId.split('_')
		return int(chromPos[1])
	
	
	# Returns whether something is in the filter or not
	def passesFilter(self, valToCheck, filterList):
		if(filterList is not None):
			if(valToCheck in filterList):
				return True
			return False
		return True
	
	
	# Returns the list of variant contexts
	def getVariantContexts(self):
		return self.variantContexts
	
	# Returns the map of variants by variant context id
	def getVariantContextsById(self):
		return self.variantContextsById
	
	# Returns the map of variant contexts ordered by sampleid
	def getVariantContextsBySample(self):
		return self.variantContextsBySample
	
	# Returns the variant contexts for a specified sampleid  or None if the sampleid is not in the map
	def getSampleVariantContexts(self, sampleid):
		if(sampleid in self.variantContextsBySample):
			return self.variantContextsBySample[sampleid]
		return None
	
	
	# Returns the sample ID of a specified variant context
	def getVariantContextSample(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getVariantContextSample()
		return None
	
	# Returns the variant context chrom of a specified variant context
	def getVariantContextChrom(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getVariantContextChrom()
		return None
	
	# Returns the variant context origin (the variant position the vacron is based on)
	def getVariantContextOrigin(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getVariantContextOrigin()
		return None
	
	# Returns the variant context origin (the variant position the vacron is based on)
	def getVariantContextOrigin2(self, varconId):
		if(varconId in self.variantContextsById):
			return str(self.variantContextsById[varconId].getVariantContextChrom())+ ":" +str(self.variantContextsById[varconId].getVariantContextOrigin())
		return None
	
	# Returns the variant context starting position of a specified variant context
	def getVariantContextStart(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getVariantContextStart()
		return None
	
	# Returns the variant context ending position of a specified variant context
	def getVariantContextEnd(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getVariantContextEnd()
		return None
	
	# Returns the acceptor variant context length for a specified variant context
	def getAcceptorContextLength(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getAcceptorContextLength()
		return None
	
	# Returns the donor variant context length for a specified variant context
	def getDonorContextLength(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getDonorContextLength()
		return None
	
	# Returns the number of acceptor reads for a specified variant context
	def getNumberOfAcceptorReads(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getNumberOfAcceptorReads()
		return None
	
	# Returns the number of donor reads for a specified variant context
	def getNumberOfDonorReads(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getNumberOfDonorReads()
		return None
	
	# Returns the acceptor/donor read ratio for a specified variant context
	def getAcceptorDonorReadRatio(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getAcceptorDonorReadRatio()
		return None
	
	# Returns the list of read identifiers for a specified variant context
	def getAcceptorReadIds(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getAcceptorReadIds()
		return None
	
	# Returns the list of donor read identifiers for a specified variant context
	def getDonorReadIds(self, varconId):
		if(varconId in self.variantContextsById):
			return self.variantContextsById[varconId].getDonorReadIds()
		return None
	
	# Returns a list with all acceptor read identifiers of all variant contexts
	def getAllAcceptorReadIds(self):
		areadIds = []
		for varconObj in self.variantContexts:
			areadIds.extend(varconObj.getAcceptorReadIds())
		return areadIds
	
	# Returns a list with all donor read identifiers of all variant contexts
	def getAllDonorReadIds(self):
		dreadIds = []
		for varconObj in self.variantContexts:
			dreadIds.extend(varconObj.getDonorReadIds())
		return dreadIds
	
	# Returns all acceptor read identifiers per variant context identifier.
	def getAllAcceptorReadIdsByVarcon(self):
		areadIds = {}
		for varconId, varconObj in self.variantContextsById.items():
			areadIds[varconId] = varconObj.getAcceptorReadIds()
		return areadIds
	
	# Returns all donor read identifiers per variant context identifier.
	def getAllDonorReadIdsByVarcon(self):
		dreadIds = {}
		for varconId, varconObj in self.variantContextsById.items():
			dreadIds[varconId] = varconObj.getDonorReadIds()
		return dreadIds
	
	# Returns the name and location of the read varcon file.
	def getVarconFileLoc(self):
		return self.variantContextFileLocation
	
	# Returns a list of all variant context ids
	def getVariantContextIds(self):
		return list(self.variantContextsById.keys())
	
	
	# Returns the number of saved variant contexts
	def getNumberOfVariantContexts(self):
		return len(self.variantContextsById)
	
	
	# Returns the list of variant context identifiers from both variant context files.
	def getVariantContextsUnion(self, otherVarconFile):
		ownVarconIds = self.getVariantContextIds()
		otherVarconIds = otherVarconFile.getVariantContextIds()
		return list(set(ownVarconIds) | set(otherVarconIds))
	
	# Returns the list of variant context identifiers present in both variant context files
	def getVariantContextsIntersect(self, otherVarconFile):
		ownVarconIds = self.getVariantContextIds()
		otherVarconIds = otherVarconFile.getVariantContextIds()
		return list(set(ownVarconIds) & set(otherVarconIds))
	
	# Returns the list of variant context identifiers in this file but not present in the other variant context file
	def getVariantContextsDifference(self, otherVarconFile):
		ownVarconIds = self.getVariantContextIds()
		otherVarconIds = otherVarconFile.getVariantContextIds()
		return list(set(ownVarconIds) - set(otherVarconIds))
	
	# Returns the list of variant context identifiers only present in one of the variant context file but not the other
	def getVariantContextsSymmetricDifference(self, otherVarconFile):
		ownVarconIds = self.getVariantContextIds()
		otherVarconIds = otherVarconFile.getVariantContextIds()
		return list(set(ownVarconIds) ^ set(otherVarconIds))
	
	
	# Returns the map containing the field number and the name
	def getVariantContextFields(self):
		return self.varconFields
	
	
	# Adds a VarconStatsFile to the current Variant Context file.
	def attachStatsFile(self, statsFileObj):
		self.variantContextStatistics = statsFileObj
		# Procedure to add the statistics to each variant context
		
	
	# Returns the attached Varcon stats file object
	def getStatsFile(self):
		return self.variantContextStatistics
	
	
	# Writes the variant context data to an output file
	def writeVariantContextFile(self, outFileLoc, sampleFilter=None, varconFilter=None, chromFilter=None):
		try:
			with open(outFileLoc, 'w') as varconOutFile:
				if(self.variantContextFileType == 'acceptor' or self.variantContextFileType == 'donor'):
					varconOutFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tNumOfReads\tReadIds\n")
				else:
					varconOutFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tAcceptorContext\tDonorContext\tAcceptorReads\tDonorReads\tADratio\tAcceptorReadsIds\tDonorReadIds\n")
				
				for varcon in self.variantContexts:
					samplePass = self.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.passesFilter(fileLineData[0], varconFilter)
					chromPass = self.passesFilter(fileLineData[2], chromFilter)
					if(samplePass and varconPass and chromPass):
						varconOutFile.write(varcon.toString()+"\n")
		except IOError as ioe:
			self.vaseUtilLogger.warning("Could not read variant contexts to " +str(ioe.filename))
	
	
	
	#====================METHODS TO OBTAIN VARIANT CONTEXT DATA BASED ON FILTERS====================
	# Returns a list/hashmap of VariantContextObjects
	def getVariantContexts(self, asList=False, varconFilter=None, sampleFilter=None, chromFilter=None):
		if(asList):
			return [x for x in self.variantContexts if(self.passesFilter(x.getVariantContextId(), varconFilter) and self.passesFilter(x.getVariantContextSample(), sampleFilter) and self.passesFilter(x.getVariantContextChrom(), chromFilter))]
		return {k:v for k,v in self.variantContextsById if(self.passesFilter(k, varconFilter) and self.passesFilter(v.getVariantContextSample(), sampleFilter) and self.passesFilter(v.getVariantContextChrom(), chromFilter))}
	
	
	# Main method that returns whether a variant (SNP or indel).
	def variantIsInContext(self, variantType, searchChrom, searchStart, searchStop):
		if(variantType=='snp'):
			return self.snpVariantIsInContext()
		if(variantType='indel'):
			return self.indelVariantIsInContext(searchChrom, searchStart, searchStop)
	
	
	# Determines whether an SNP variant is located in an already existing variant context.
	def snpVariantIsInContext(self, varchrom, varpos):
		for varcon in self.variantContexts:
			if(varchrom == varcon.getVariantContextChrom()):
				if(vcfVarPos >= varcon.getVariantContextStart() and vcfVarPos <= varcon.getVariantContextEnd()):
					return True
		return False
	
	
	# Determines whether an indel variant is located within an existing variant context (indelLeftPos and indelRightPos can be the used search window)
	def indelVariantIsInContext(self, indelChrom, indelLeftPos, indelRightPos):
		for varcon in self.variantContexts:
			if(indelChrom == varcon.getVariantContextChrom()):
				if(indelLeftPos <= varcon.getContextStart() and indelRightPos >= varcon.getContextStart()):
					return True
				if(indelLeftPos <= varcon.getContextEnd() and indelRightPos >= varcon.getContextEnd()):
					return True
				if(indelLeftPos >= varcon.getContextStart() and indelRightPos <= varcon.getContextEnd()):
					return True
		return False
	
	
	# Adds the variant context object
	def addVariantContext(self, varconType, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconReads, otherVarconReads=None, varconALength=None, varconDLength=None):
		varconObj = VariantContext(varconType, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconReads, otherVarconReads, varconALength, varconDLength)
		self.variantContexts.append(varconObj)
		self.variantContextsById[varconId] = varconObj
		if(varconSample not in self.variantContextsBySample):
			self.variantContextsBySample[varconSample] = []
		self.variantContextsBySample[varconSample].append(varconObj)
