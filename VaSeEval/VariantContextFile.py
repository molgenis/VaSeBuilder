#!/usr/bin/env python
from VariantContext import VariantContext
import logging

class VariantContextFile:
	# Constructor
	def __init__(self, varconFileLoc):
		self.variantContexts = {}	# Will save VariantContext objects for all variants.
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
		self.readVariantContextFile(varconFileLoc)	# Read the provided variant context file
	
	
	# Reads the variant context file.
	def readVariantContextFile(self, fileLoc):
		try:
			with open(fileLoc, 'r') as vcFile:
				for fileLine in vcFile.readlines():
					fileLine = fileLine.strip()
					
					# Process the variant context data. The check ensures the header line is not processed.
					if(not fileLine.startswith("Variant")):
						fileLineData = fileLine.split("\t")
						self.variantContexts[fileLineData[0]] = VariantContext(fileLineData[0], fileLineData[1], fileLineData[2], fileLineData[3])
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not read variant context file " +str(fileLoc))
	
	
	# Returns a specified variant context. Returns None if the variant is not available.
	def getVariantContext(self, variantId):
		if(self.hasVariantContext(variantId)):
			return self.variantContexts[variantId]
		else:
			self.vaseEvalLogger.debug("Variant " +str(variantId)+ " was not found.")
			return None
	
	
	#Returns whether a variant is in the variant context file.
	def hasVariantContext(self, variantId):
		if(variantId in self.variantContexts):
			return True
		return False
	
	
	# Returns all variant contexts.
	def getVariantContexts(self):
		return self.variantContexts
	
	
	# Adds a donor read identifier to the specified variant context.
	def addDonorReadToVariantContext(self, variantId, donorReadId):
		if(variantId in self.variantContexts):
			self.variantContexts[variantId].addDonorReadIdentifier(donorReadId)
			return True
		else:
			self.vaseEvalLogger.debug("Context for variant " +str(variantId)+ " was not found and the donor read identifier could therefore not be added.")
			return False
	
	
	# Adds the list of donor read identifiers to the variant context.
	def addDonorReadsToVariantContext(self, variantId, donorReadIds):
		if(variantId in self.variantContexts):
			self.variantContexts[variantId].addDonorReadIdentifiers(donorReadIds)
			return True
		else:
			self.vaseEvalLogger.debug("Context for variant " +str(variantId)+ " was not found and could therefore not add donor read identifiers.")
			return False
	
	
	# Adds a donor read identifier to the specified variant context.
	def addAcceptorReadToVariantContext(self, variantId, acceptorReadId):
		if(variantId in self.variantContexts):
			self.variantContexts[variantId].addAcceptorReadIdentifier(acceptorReadId)
			return True
		else:
			self.vaseEvalLogger.debug("Context for variant " +str(variantId)+ " was not found and the acceptor read identifier could therefore not be added.")
			return False
	
	
	# Adds the list of acceptor read identifiers to the variant context.
	def addAcceptorReadsToVariantContext(self, variantId, acceptorReadIds):
		if(variantId in self.variantContexts):
			self.variantContexts[variantId].addAcceptorReadIdentifiers(acceptorReadIds)
			return True
		else:
			self.vaseEvalLogger.debug()
			return False
	
	
	# Returns the number of contexts in the variant context file.
	def getNumberOfContexts(self):
		return len(self.variantContexts)
