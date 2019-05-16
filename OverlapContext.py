import statistics
from DonorBamRead import DonorBamRead


class OverlapContext:
    # Saves the data associated with the context overlapping with the
    # variant.
    def __init__(self, variantid, sampleid, ovconchrom, ovconorigin,
                 ovconstart, ovconend, bamreads):
        self.contextId = variantid
        self.sampleId = sampleid
        self.contextChrom = ovconchrom
        self.contextOrigin = ovconorigin
        self.contextStart = ovconstart
        self.contextEnd = ovconend
        self.contextBamReads = bamreads
        self.unmappedReadMateIds = []

    # ===METHODS TO GET DATA OF THE OVERLAP CONTEXT============================
    # Returns the context identifier.
    def getContextId(self):
        return self.contextId

    # Returns the sample identifier the context was based on.
    def getSampleId(self):
        return self.sampleId

    # Returns the context chromosome name the context is located on.
    def getContextChrom(self):
        return self.contextChrom

    # Returns the origin position of the context.
    def getContextOrigin(self):
        return self.contextOrigin

    # Returns the context start position.
    def getContextStart(self):
        return self.contextStart

    # Returns the end position of the context.
    def getContextEnd(self):
        return self.contextEnd

    # Returns the bam reads associated with the context as a list of
    # BamRead objects.
    def getContextBamReads(self):
        return self.contextBamReads

    # Returns the list of BAM read ids that have an unmapped mate.
    def getUnmappedReadMateIds(self):
        return self.unmappedReadMateIds

    # ===METHODS TO GET DATA (REQUIRING SOME CALCULATION) FROM THE=============
    # ===OVERLAP CONTEXT===
    # Returns the lengt of the context.
    def getContextLength(self):
        return abs(self.contextEnd - self.contextStart)

    # Returns the distance of the context start from the context origin.
    def getStartDistanceFromOrigin(self):
        return abs(self.contextOrigin - self.contextStart)

    # Returns the distance of the context end from the context origin.
    def getEndDistanceFromOrigin(self):
        return abs(self.contextEnd - self.contextOrigin)

    # ===METHODS TO OBTAIN CONTEXT READ INFORMATION============================
    # Returns the number of saved context reads.
    def getNumberOfContextReads(self):
        return len(self.contextBamReads)

    # Returns a list of BAM read identifiers in the current context.
    def getContextBamReadIds(self):
        return [x.getBamReadId() for x in self.contextBamReads]

    # Returns a list of all left positions for all BAM reads.
    def getContextBamReadStarts(self):
        return [x.getBamReadRefPos() for x in self.contextBamReads]

    # Returns a list of all left positions for all R1 BAM reads.
    def getContextBamReadLeftPositions(self):
        return [x.getBamReadRefPos()
                for x in self.contextBamReads if (x.isRead1())]

    # Returns a list of BAM read ending positions for all BAM reads.
    def getContextBamReadEnds(self):
        return [x.getBamReadRefEnd() for x in self.contextBamReads]

    # Returns a list of all right positions for all R2 BAM reads.
    def getContextBamReadRightPositions(self):
        return [x.getBamReadRefEnd()
                for x in self.contextBamReads if (x.isRead2())]

    # Returns a list of all lengths for all BAM reads.
    def getContextBamReadLengths(self):
        return [x.getBamReadLength() for x in self.contextBamReads]

    # Returns a list of BAM read sequences in the current context.
    def getContextBamReadSeqs(self):
        return [x.getBamReadSequence() for x in self.contextBamReads]

    # Returns a list of qualities of all BAM reads.
    def getContextBamReadQualities(self):
        return [x.getBamReadQual() for x in self.contextBamReads]

    # Returns a list of Q-scores of all BAM reads.
    def getContextBamReadQScores(self):
        return [x.getBamReadQScores() for x in self.contextBamReads]

    # Returns a list of all BAM read MapQ values.
    def getContextBamReadMapQs(self):
        return [x.getMappingQual() for x in self.contextBamReads]

    # Returns whether a BAM read is in the context based on the provided
    # read identifier.
    def readIsInContext(self, readId):
        return readId in self.getContextBamReadIds()

    # ===METHODS TO ADD/SET CONTEXT DATA=======================================
    # Adds the read id of a BAM read with an unmapped mate.
    def addUnmappedMateId(self, uReadId):
        self.unmappedReadMateIds.append(uReadId)

    # Sets the list of unmapped read mate ids.
    def setUnmappedMateIds(self, mateIds):
        self.unmappedReadMateIds = mateIds

    # Returns whether a BAM read in the context has an unmapped mate.
    def readHasUnmappedMate(self, readId):
        return readId in self.unmappedReadMateIds

    # ===STATISTICS METHODS FOR A VARIANT CONTEXT==============================
    # Returns the average and median read length.
    def getAverageAndMedianReadLength(self):
        avgMedLen = []
        for contextread in self.contextBamReads:
            avgMedLen.append(contextread.getBamReadLength())
        return ([statistics.mean(avgMedLen), statistics.median(avgMedLen)])

    # Returns the average and median read quality.
    def getAverageAndMedianReadQual(self):
        avgMedQual = []
        for contextread in self.contextBamReads:
            avgMedQual.append(contextread.getAverageQscore())
        return ([statistics.mean(avgMedQual), statistics.median(avgMedQual)])

    # Returns the average and median read MapQ of this variant context.
    def getAverageAndMedianReadMapQ(self):
        avgMedMapQ = []
        for contextread in self.contextBamReads:
            avgMedMapQ.append(contextread.getMappingQual())
        return ([statistics.mean(avgMedMapQ), statistics.median(avgMedMapQ)])

    # ===SOME OTHER METHODS====================================================
    # Returns a string representation of the overlap context.
    def toString(self):
        return (str(self.contextId) + "\t"
                + str(self.sampleId) + "\t"
                + str(self.contextChrom) + "\t"
                + str(self.contextOrigin) + "\t"
                + str(self.contextStart) + "\t"
                + str(self.contextEnd) + "\t"
                + str(len(self.contextBamReads)) + "\t"
                + ';'.join([x.getBamReadId() for x in self.contextBamReads]))

    # Returns a statistics string representation of the overlap context.
    def toStatisticsString(self):
        avgMedLens = self.getAverageAndMedianReadLength()
        avgMedQuals = self.getAverageAndMedianReadQual()
        avgMedMapQ = self.getAverageAndMedianReadMapQ()
        return (f"{self.contextId}\t{avgMedLens[0]}\t{avgMedLens[1]}\t"
                f"{avgMedQuals[0]}\t{avgMedQuals[1]}\t{avgMedMapQ[0]}\t"
                f"{avgMedMapQ[1]}")

    # Compares the current OverlapContext to another OverlapContext and
    # returns the differences.
    def compare(self, otherOverlapContext):
        differences = {}
        if (self.contextId != otherOverlapContext.getContextId()):
            differences[1] = [self.contextId,
                              otherOverlapContext.getContextId()]
        if (self.sampleId != otherOverlapContext.getSampleId()):
            differences[2] = [self.sampleId,
                              otherOverlapContext.getSampleId()]
        if (self.contextChrom != otherOverlapContext.getContextChrom()):
            differences[3] = [self.contextChrom,
                              otherOverlapContext.getContextChrom()]
        if (self.contextOrigin != otherOverlapContext.getContextOrigin()):
            differences[4] = [self.contextOrigin,
                              otherOverlapContext.getContextOrigin()]
        if (self.contextStart != otherOverlapContext.getContextStart()):
            differences[5] = [self.contextStart,
                              otherOverlapContext.getContextStart()]
        if (self.contextEnd != otherOverlapContext.getContextEnd()):
            differences[6] = [self.contextEnd,
                              otherOverlapContext.getContextEnd()]
        if (self.getContextBamReadIds().sort() != otherOverlapContext.getContextBamReadIds().sort()):
            differences[7] = [self.contextBamReads,
                              otherOverlapContext.getContextBamReads()]
        return differences
