import logging
from VariantContext import VariantContext

class VariantContextFile:
	def __init__(self, vuhelper, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = vuhelper
		self.variantContextsBySample = {}
		self.variantContextsById = {}
		self.variantContexts = []
		self.readVariantContextFile(self, fileLoc, sampleFilter, varconFilter, chromFilter, posFilter)
		self.varconFields = {1: 'variant context id',
			2 : 'sample id',
			3 : 'chromosome',
			4 : 'origin',
			5 : 'start pos',
			6 : 'end pos'
			}
	
	
	# Reads the varcon files and saves data according to set filters
	def readVariantContextFile(self, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		try:
			with open(fileLoc, 'r') as vcFile:
				next(vcFile)	# Skip the header line
				for fileLine in vcFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.vuh.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.vuh.passesFilter(fileLineData[0], varconFilter)
					chromPass = self.vuh.passesFilter(fileLineData[2], chromFilter)
					posPass = self.vuh.passesFilter(self.getVariantContextVarPos(fileLineData[0]), posFilter)
					
					if(samplePass and chromPass and posPass):
						varconObj = VariantContext(fileLineData[0], fileLineData[1], fileLineData[2], int(fileLineData[3]), int(fileLineData[4]), int(fileLineData[5]))
						if(fileLineData[1] not in self.variantContextsBySample):
							self.variantContextsBySample[fileLineData[1]] = []
						self.variantContextsBySample[fileLineData[1]].append(varconObj)
						self.variantContextsById[fileLineData[0]] = varconObj
						self.variantContexts.append(varconObj)
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read varcon file")
			exit()
	
	
	# Compare the variant contexts to another set
	def compare(self, varconFile):
		varconNotInOther = 0
		otherVarcon = varconFile.getVariantContextsById()
		
		# Start the comparison.
		for varconId, varconObj in self.variantContextsById.items():
			if(varconId in otherVarcon):
				vcdiffs = varconObj.compare(otherVarcon[varconId])
				if(len(vcdiffs)>0):
					vulmsg = "Variant context " +str(varconId)+ " differs on "
					#self.vaseUtilLogger.info("Variant context " +str(varconId)+ " differs")
					
					diffields = []
					for vd in vcdiffs:
						if(vd in self.varconFields):
							diffields.append(self.varconFields[vd])
					vulmsg += ", ".join(diffields)
					
					print(vulmsg)
			else:
				self.vaseUtilLogger.info("Variant Context " +str(varconId)+ "not found in other")
	
	
	# Returns the variant position based on the SNP##_#### identifier
	def getVariantContextVarPos(self):
		if(self.variantContextId.startswith('SNP')):
			chromPos = self.variantContextId.split('_')
			return int(chromPos[1])
		return -1
	
	
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
