import statistics
from OverlapContext import OverlapContext
from DonorBamRead import DonorBamRead

class VariantContext:
	# Sets the variant context data.
	def __init__(self, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, acceptorReads, donorReads, acceptorContext=None, donorContext=None):
		self.contextId = varconId
		self.sampleId = varconSample
		self.variantContextChrom = varconChrom
		self.variantContextOrigin = varconOrigin
		self.variantContextStart = varconStart
		self.variantContextEnd = varconEnd
		self.variantContextAReads = acceptorReads	# Saves the acceptor reads that overlap with the entire variant context
		self.variantContextDReads = donorReads	# Saves the donor reads that overlap with the entire variant context
		self.variantAcceptorContext = acceptorContext	# Saves the acceptor context (this is the context from acceptor reads overlapping with the variant itself)
		self.variantDonorContext = donorContext	# Saves the donor context (this is the context from donor reads overlapping with the variant itself)
		self.unmappedAcceptorMateIds = []
		self.unmappedDonorMateIds = []
	
	
	
	# ====================METHODS TO OBTAIN DATA FROM THE VARIANT CONTEXT DATA====================
	# Returns the variant context identifier
	def getVariantContextId(self):
		return self.contextId
	
	# Returns the variant context sample
	def getVariantContextSample(self):
		return self.sampleId
	
	# Returns the variant context chromosome
	def getVariantContextChrom(self):
		return self.variantContextChrom
	
	# Returns the variant context origin
	def getVariantContextOrigin(self):
		return self.variantContextOrigin
	
	# Returns the variant context start position
	def getVariantContextStart(self):
		return self.variantContextStart
	
	# Returns the variant context end position
	def getVariantContextEnd(self):
		return self.variantContextEnd
	
	# Returns the variant context acceptor reads
	def getAcceptorReads(self):
		return self.variantContextAReads
	
	# Returns the variant context donor reads
	def getDonorReads(self):
		return self.variantContextDReads
	
	# Returns the acceptor context (the context from acceptor reads overlapping the variant itself)
	def getAcceptorContext(self):
		return self.variantAcceptorContext
	
	# Returns the donor context (the context from donor reads overlapping the variant itself)
	def getDonorContext(self):
		return self.variantDonorContext
	
	# Returns a list of acceptor reads overlapping with the variant context
	def getUnmappedAcceptorMateIds(self):
		return self.unmappedAcceptorMateIds
	
	# Returns a list of donor reads overlapping with the variant context
	def getUnmappedDonorMateIds(self):
		return self.unmappedDonorMateIds
	
	
	
	# ====================METHODS TO GET DATA (REQUIRING SOME CALCULATING) OF THE VARIANT CONTEXT====================
	# Returns the variant context length
	def getVariantContextLength(self):
		return abs(self.variantContextEnd - self.variantContextStart)
	
	# Returns the distance of the variant context start position from the variant context origin
	def getStartDistanceFromOrigin(self):
		return abs(self.variantContextOrigin - self.variantContextStart)
	
	# Returns the distance of the variant context end position from the variant context origin
	def getEndDistanceFromOrigin(self):
		return abs(self.variantContextEnd - self.variantContextOrigin)
	
	
	
	# ====================METHODS TO OBTAIN VARIANT CONTEXT ACCEPTOR READ DATA====================
	# Returns the number of variant context acceptor reads
	def getNumberOfAcceptorReads(self):
		return len(self.variantContextAReads)
	
	# Returns the identifiers of acceptor reads overlapping with the variant context
	def getAcceptorReadIds(self):
		return [x.getBamReadId() for x in self.variantContextAReads]
	
	# Returns the list of left most acceptor read positions
	def getAcceptorReadStarts(self):
		return [x.getBamReadRefPos() for x in self.variantContextAReads]
	
	# Returns a list of the left most positions of the R1 variant context accecptor BAM reads
	def getAcceptorReadLeftPositions(self):
		return [x.getBamReadRefPos() for x in self.variantContextAReads if x.isRead1()]
	
	# Returns the list of all end positions for all variant context acceptor reads 
	def getAcceptorReadEnds(self):
		return [x.getBamReadRefEnd() for x in self.variantContextAReads]
	
	# Returns a list of the right most positions ofr the R2 variant context acceptor BAM reads
	def getAcceptorReadRightPositions(self):
		return [x.getBamReadRefEnd() for x in self.variantContextAReads if x.isRead2()]
	
	
	
	# ====================METHODS TO OBTAIN VARIANT CONTEXT DONOR READ DATA====================
	# Returns the number of variant context donor reads
	def getNumberOfDonorReads(self):
		return len(self.variantContextDReads)
	
	# Returns the identifiers of donor reads overlapping with the variant context
	def getDonorReadIds(self):
		return [x.getBamReadId() for x in self.variantContextDReads]
	
	# Returns the list of variant context donor read starting positions
	def getDonorReadStarts(self):
		return [x.getBamReadRefPos() for x in self.variantContextDReads]
	
	# Returns the list of variant context donor read pairs left most positions (start pos of read 1)
	def getDonorReadLeftPositions(self):
		return [x.getBamReadRefPos() for x in self.variantContextDReads if(x.isRead1())]
	
	# Returns a list of all donor read ending positions
	def getDonorReadEnds(self):
		return [x.getBamReadRefEnd() for x in self.variantContextDReads]
	
	# Returns a list of all variant context donor reads right most positions (end pos of read 2)
	def getDonorReadRightPositions(self):
		return [x.getBamReadRefEnd() for x in self.variantContextDReads if(x.isRead2())]
	
	
	
	# ====================METHODS TO ADD DATA TO THE VARIANT CONTEXT====================
	# Sets the acceptor context of the variant context
	def setAcceptorContext(self, acceptorContext):
		self.variantAcceptorContext = acceptorContext
	
	# Sets the donor context of the variant context
	def setDonorContext(self, donorContext):
		self.variantDonorContext = donorContext
	
	# Adds an acceptor context to the variant context by creating it from data values
	def addAcceptorContext(self, contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, acceptorReads):
		self.variantAcceptorContext = OverlapContext(contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, acceptorReads)
	
	# Adds a donor context to the variant context by creating it from data values
	def addDonorContext(self, contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, donorReads):
		self.variantDonorContext = OverlapContext(contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, donorReads)
	
	
	
	# ====================METHODS TO OBTAIN VARIANT CONTEXT UNMAPPED MATE READ DATA====================	
	# Returns the variant context acceptor read ids that have an unmapped mate
	def getUnmappedAcceptorReadIds(self):
		return unmappedAcceptorMateIds
	
	# Returns the variant context donor read ids that have an unmapped mate
	def getUnmappedDonorReadIds(self):
		return unmappedDonorMateIds
	
	# Adds a variant context appector mate identifier
	def addUnmappedAcceptorMateId(sef, mateId):
		self.unmappedAcceptorMateIds.append(mateId)
	
	# Adds a variant context donor mate identifier
	def addUnmappedDonorMateId(self, mateId):
		self.unmappedDonorMateIds.append(mateId)
	
	# Sets the variant context unmapped acceptor mate ids
	def setUnmappedAcceptorMateIds(self, mateIds):
		self.unmappedAcceptorMateIds = mateIds
	
	# Sets the variant context unmapped donor mate ids
	def setUnmappedDonorMateIds(self, mateIds):
		self.unmappedDonorMateIds = mateIds
	
	# Returns whether a specified variant context acceptor read has an unmapped mate
	def acceptorReadHasUnmappedMate(self, readId):
		return readId in self.unmappedAcceptorMateIds
	
	# Returns whether a specified variant context donor read has an unmapped mate
	def donorReadHasUnmappedMate(self, readId):
		return readId in self.unmappedDonorMateIds
	
	# Returns the number of variant context acceptor reads with an unmapped mate
	def getNumberOfUnmappedAcceptorMates(self):
		return len(self.unmappedAcceptorMateIds)
	
	# Returns the number of variant context donor reads with an unmapped mate
	def getNumberOfUnmappedDonorMates(self):
		return len(self.unmappedDonorMateIds)
	
	
	
	# ====================METHODS TO ADD UNMAPPED MATES TO THE ACCEPTOR AND DONOR CONTEXT====================
	# Sets the unmapped mate ids for the acceptor context
	def setAcceptorContextUnmappedMates(self, mateIds):
		self.variantAcceptorContext.setUnmappedMateIds(mateIds)
	
	# Adds an unmapped read id to the acceptor context
	def addAcceptorContextUnmappedMate(self, uReadId):
		self.variantAcceptorContext.addUnmappedMateId(ureadId)
	
	# Sets the unmapped mate ids for the donor context
	def setDonorContextUnmappedMates(self, mateIds):
		self.variantDonorContext.setUnmappedMateIds(mateIds)
	
	# Adds an unmapped read id to the donor context
	def addDonorContextUnmappedMateId(self, uReadId):
		self.variantDonorContext.setUnmappedMateId(uReadId)
	
	
	
	# ====================METHODS TO OBTAIN STATISTICS OF THE VARIANT CONTEXT====================
	# Returns the average and median quality of the acceptor reads associated with
	def getAverageAndMedianAcceptorReadQual(self):
		return self.getAverageAndMedianReadQual(self.variantContextAReads)
	
	# Returns the average and median quality of the donor reads associated with this variant context
	def getAverageAndMedianDonorReadQual(self):
		return self.getAverageAndMedianReadQual(self.variantContextDReads)
	
	# Returns the average and median read quality 
	def getAverageAndMedianReadQual(self, contextReads):
		if(contextReads is not None):
			avgMedQual = []
			for contextread in contextReads:
				avgMedQual.append(contextread.getAverageQscore())
			return ([statistics.mean(avgMedQual), statistics.median(avgMedQual)])
		return None
	
	
	# Returns the average and median mapq values of the acceptor reads associated with this variant context
	def getAverageAndMedianAcceptorReadMapQ(self):
		return self.getAverageAndMedianReadMapQ(self.variantContextAReads)
	
	# Returns the average and median mapq values of the donor reads associated with this variant context
	def getAverageAndMedianDonorReadMapQ(self):
		return self.getAverageAndMedianReadMapQ(self.variantContextDReads)
	
	# Returns the average and median read MapQ of this variant context
	def getAverageAndMedianReadMapQ(self, contextReads):
		if(contextReads is not None):
			avgMedMapQ = []
			for contextread in contextReads:
				avgMedMapQ.append(contextread.getMappingQual())
			return ([statistics.mean(avgMedMapQ), statistics.median(avgMedMapQ)])
		return None
	
	
	# Returns the average and median length of the acceptor reads associated with this variant context
	def getAverageAndMedianAcceptorReadLength(self):
		return self.getAverageAndMedianReadLength(self.variantContextAReads)
	
	# Returns the average and median length of the donor reads associated with this variant context
	def getAverageAndMedianDonorReadLength(self):
		return self.getAverageAndMedianReadLength(self.variantContextDReads)
	
	# Returns the average and median read length
	def getAverageAndMedianReadLength(self, contextReads):
		if(contextReads is not None):
			avgMedLen = []
			for contextread in contextReads:
				avgMedLen.append(contextread.getBamReadLength())
			return ([statistics.mean(avgMedLen), statistics.median(avgMedLen)])
		return None
	
	
	
	
	
	# ====================METHODS TO OBTAIN ACCEPTOR CONTEXT DATA====================
	# Returns the acceptor context identifier (should be the same as the variant context id)
	def getAcceptorContextId(self):
		return self.variantAcceptorContext.getContextId()
	
	# Returns the acceptor context sample id (should be the same as the variant context sample id)
	def getAcceptorContextSampleId(self):
		return self.variantAcceptorContext.getSampleId()
	
	# Returns the chromosome of the acceptor context
	def getAcceptorContextChrom(self):
		return self.variantAcceptorContext.getContextChrom()
	
	# Returns the origin position of the acceptor context
	def getAcceptorContextOrigin(self):
		return self.variantAcceptorContext.getContextOrigin()
	
	# Returns the starting position of the acceptor context
	def getAcceptorContextStart(self):
		return self.variantAcceptorContext.getContextStart()
	
	# Returns the ending position of the acceptor context
	def getAcceptorContextEnd(self):
		return self.variantAcceptorContext.getContextEnd()
	
	# Returns the length of the acceptor context
	def getAcceptorContextLength(self):
		return self.variantAcceptorContext.getContextLength()
	
	# Returns the acceptor context reads
	def getAcceptorContextReads(self):
		return self.variantAcceptorContext.getContextBamReads()
	
	# Returns the lost of read ids in the acceptor context
	def getAcceptorContextReadIds(self):
		return self.variantAcceptorContext.getContextBamReadIds()
	
	# Returns a list of all acceptor context BAM read start positions
	def getAcceptorContextReadStarts(self):
		return self.variantAcceptorContext.getContextBamReadStarts()
	
	# Returns a list of all acceptor context R1 BAM read left positons
	def getAcceptorContextReadLeftPositions(self):
		return self.variantAcceptorContext.getContextBamReadLeftPositions()
	
	# Returns a list of all acceptor context BAM read end positions
	def getAcceptorContextReadEnds(self):
		return self.variantAcceptorContext.getContextBamReadEnds()
	
	# Returns a list of all acceptor context R2 BAM read end positions
	def getAcceptorContextReadRightPositions(self):
		return self.variantAcceptorContext.getContextBamReadRightPositions()
	
	# Returns a list of all acceptor context BAM read lengths
	def getAcceptorContextReadLengths(self):
		return self.variantAcceptorContext.getContextBamReadLengths()
	
	# Returns the list of acceptor context unmapped mate read ids
	def getAcceptContextUnmappedMateIds(self):
		return self.variantAcceptorContext.getUnmappedReadMateIds()
	
	
	
	# ====================METHODS TO OBTAIN DONOR CONTEXT DATA====================
	# Returns the donor context identifier (should be the same as the variant context id)
	def getDonorContextId(self):
		return self.variantDonorContext.getContextId()
	
	# Returns the acceptor context sample id (should be the same as the variant context sample id)
	def getDonorContextSampleId(self):
		return self.variantDonorContext.getSampleId()
	
	# Returns the chromosome of the donor context
	def getDonorContextChrom(self):
		return self.variantDonorContext.getContextChrom()
	
	# Returns the origin position of the donor context
	def getDonorContextOrigin(self):
		return self.variantDonorContext.getContextOrigin()
	
	# Returns the starting position of the donor context
	def getDonorContextStart(self):
		return self.variantDonorContext.getContextStart()
	
	# Returns the ending position of the donor context
	def getDonorContextEnd(self):
		return self.variantDonorContext.getContextEnd()
	
	# Returns the length of the acceptor context
	def getDonorContextLength(self):
		return self.variantDonorContext.getContextLength()
	
	# Returns the donor context reads
	def getDonorContextReads(self):
		return self.variantDonorContext.getContextBamReads()
	
	# Returns the donor context read identifers
	def getDonorContextReadIds(self):
		return self.variantDonorContext.getContextBamReadIds()
	
	# Returns a list of donor context read start positions
	def getDonorContextReadStarts(self):
		return self.variantDonorContext.getContextBamReadStarts()
	
	# Returns a list of all acceptor context R1 BAM read left positons
	def getDonorContextReadLeftPositions(self):
		return self.variantDonorContext.getContextBamReadLeftPositions()
	
	# Returns a list of donor context read end positions
	def getDonorContextReadEnds(self):
		return self.variantDonorContext.getContextBamReadEnds()
	
	# Returns a list of all acceptor context R2 BAM read right positons
	def getDonorContextReadRightPositions(self):
		return self.variantDonorContext.getContextBamReadRightPositions()
	
	# Returns a list of donor context BAM read lengths
	def getDonorContextReadLengths(self):
		return self.variantDonorContext.getContextBamReadLengths()
	
	# Returns the list of donor context unmapped mate read ids
	def getDonorContextUnmappedMateIds(self):
		return self.variantDonorContext.getUnmappedReadMateIds
	
	
	
	# ====================METHODS TO PRODUCE SOME OUTPUT ABOUT THE VARIANT CONTEXT====================
	# Returns a varcon.txt string representation of the variant context
	def toString(self):
		return str(self.contextId)+ "\t" +str(self.sampleId)+ "\t" +str(self.variantContextChrom)+ "\t" +str(self.variantContextOrigin)+ "\t" +str(self.variantContextStart)+ "\t" +str(self.variantContextEnd)+ "\t"  +str(self.variantAcceptorContext.getContextLength())+ "\t" +str(self.variantDonorContext.getContextLength())+ "\t" +str(len(self.variantContextAReads))+ "\t" +str(len(self.variantContextDReads))+ "\t"  +str(float(len(self.variantContextAReads)/len(self.variantContextDReads)))+ "\t" +';'.join(self.getAcceptorReadIds())+ "\t" +';'.join(self.getDonorReadIds())
	
	# Returns a varconstats.txt string representation of the variant context
	def toStatisticsString(self):
		aReadLen = self.getAverageAndMedianAcceptorReadLength()
		dReadLen = self.getAverageAndMedianDonorReadLength()
		aReadQual = self.getAverageAndMedianAcceptorReadQual()
		dReadQual = self.getAverageAndMedianDonorReadQual()
		aReadMapQ = self.getAverageAndMedianAcceptorReadMapQ()
		dReadMapQ = self.getAverageAndMedianDonorReadMapQ()
		return str(self.contextId)+ "\t" +str(aReadLen[0])+ "\t" +str(dReadLen[0])+ "" +str(aReadLen[1])+ "" +str(dReadLen[1])+ "\t" +str(aReadQual[0])+ "\t" +str(dReadQual[0])+ "\t" +str(aReadQual[1])+ "\t" +str(dReadQual[1])+ "\t" +str(aReadMapQ[0])+ "\t" +str(dReadMapQ[0])+ "\t" +str(aReadMapQ[1])+ "\t" +str(dReadMapQ[1])
