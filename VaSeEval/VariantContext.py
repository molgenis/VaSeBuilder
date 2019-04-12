#!/usr/bin/env python
import logging

class VariantContext:
	# Constructor saving the variant identifier and context chromosome, starting and ending positions.
	def __init__(self, varId, conChr, conStart, conEnd):
		self.variantId = varId
		self.contextChrom = conChr
		self.contextStart = conStart
		self.contextEnd = conEnd
		self.contextDonorReadIds = []
		self.contextAcceptorReadIds = []
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
	
	
	# Returns the identifier of the variant.
	def getVariantId(self):
		return self.variantId
	
	
	# Returns the chromosome of the variant context.
	def getContextChrom(self):
		return self.contextChrom
	
	
	# Returns the starting position of the variant context.
	def getContextStart(self):
		return self.contextStart
	
	
	# Returns the ending position of the variant context.
	def getContextEnd(self):
		return self.contextEnd
	
	
	# Returns the length of the context.
	def getContextLength(self):
		return (self.contextEnd - self.contextStart)
	
	
	# Returns the donor read identifiers associated with the variant context.
	def getContextDonorReadIds(self):
		return self.contextDonorReadIds
	
	
	# Returns the acceptor read identifiers associated with the variant context.
	def getContextAcceptorReadIds(self):
		return self.contextAcceptorReadIds
	
	
	# Returns a donor read id based on read id. 
	def getContextDonorRead(self, readId):
		try:
			return self.contextAcceptorReadIds[self.contextAcceptorReadIds.index(readId)]
		except ValueError as vae:
			self.vaseEvalLogger.debug("Could not find donor read " +str(readId)+ " for context of variant " +str(self.variantId))
			return None
	
	
	# Returns an acceptor read id based on read id.
	def getContextAcceptorRead(self, readId):
		try:
			return self.contextDonorReadIds[self.contextDonorReadIds.index(readId)]
		except ValueError as vae:
			self.vaseEvalLogger.debug("Could not find acceptor read " +str(readId)+ " for context of variant " +str(self.variantId))
			return None
	
	
	# Returns the number of donor reads associated with the variant context.
	def getNumberOfDonorReads(self):
		return len(self.contextDonorReadIds)
	
	
	# Returns the number of acceptor reads associated with the variant context.
	def getNumberOfAcceptorReads(self):
		return len(self.contextAcceptorReadIds)
	
	
	# Returns whether a specified donor read is located in the variant context.
	def donorReadIsInVariantContext(self, readId):
		return readId in self.contextDonorReadIds
	
	
	# Returns whether a specified acceptor read is located in the variant context.
	def acceptorReadIsInVariantContext(self, readId):
		return readId in self.contextAcceptorReadIds
	
	
	# Adds a read identifier to the list of donor read identifiers associated with the variant context.
	def addDonorReadIdentifier(self, readId):
		if(type(readId) is not list):
			self.contextDonorReadIds.append(readId)
		else:
			self.addDonorReadIdentifiers(readId)
	
	
	# Adds a list of donor read identifiers to the list of donor read identifiers associated with the variant context.
	def addDonorReadIdentifiers(self, listOfReadIds):
		self.contextDonorReadIds.extend(listOfReadIds)
	
	
	# Adds a read identifier to the list of acceptor read identifiers associated with the variant context.
	def addAcceptorReadIdentifier(self, readId):
		if(type(readId) is not list):
			self.contextAcceptorReadIds.append(readId)
		else:
			self.addAcceptorReadIdentifiers(readId)
	
	
	# Adds a list of donor read identifiers to the list of acceptor read identifiers associated with the variant context.
	def addAcceptorReadIdentifiers(self, listOfReadIds):
		self.contextAcceptorReadIds.extend(listOfReadIds)
