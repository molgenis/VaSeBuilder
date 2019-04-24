import logging
from VariantContext import VariantContext

class VariantContextFile:
	def __init__(self, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.variantContextFileLocation = fileLoc
		self.variantContextsBySample = {}
		self.variantContextsById = {}
		self.variantContexts = []
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
		self.readVariantContextFile(self, fileLoc, sampleFilter, varconFilter, chromFilter, posFilter)	# Read the provided variant context file with set optional filters
	
	
	# Reads the varcon files and saves data according to set filters
	def readVariantContextFile(self, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		try:
			with open(fileLoc, 'r') as vcFile:
				next(vcFile)	# Skip the header line
				for fileLine in vcFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.passesFilter(fileLineData[0], varconFilter)
					chromPass = self.passesFilter(fileLineData[2], chromFilter)
					posPass = self.passesFilter(self.getVariantContextVarPos(fileLineData[0]), posFilter)
					
					if(samplePass and chromPass and posPass):
						varconObj = VariantContext(fileLineData[0], fileLineData[1], fileLineData[2], int(fileLineData[3]), int(fileLineData[4]), int(fileLineData[5]), int(fileLineData[6]), int(fileLineData[7]), int(fileLineData[8]), int(fileLineData[9]), float(fileLineData[10]), fileLineData[11].split(';'), fileLineData[12].split(';'))
						if(fileLineData[1] not in self.variantContextsBySample):
							self.variantContextsBySample[fileLineData[1]] = []
						self.variantContextsBySample[fileLineData[1]].append(varconObj)
						self.variantContextsById[fileLineData[0]] = varconObj
						self.variantContexts.append(varconObj)
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read varcon file")
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
	def getNumberOfDonorReads(self):
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
