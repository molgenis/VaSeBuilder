import logging
import statistics
from OverlapContext import OverlapContext
from VariantContext import VariantContext

class VariantContextFile:
    def __init__(self, fileLoc=None, sampleFilter=None, varconFilter=None, chromFilter=None):
        self.vaseUtilLogger = logging.getLogger("VaSe_Logger")
        self.variantContextFileLocation = fileLoc
        self.variantContexts = {}
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
            self.readVariantContextFile(self, fileLoc, sampleFilter, varconFilter, chromFilter, posFilter)    # Read the provided variant context file with set optional filters



    # ====================METHODS TO GET DATA FROM THE VARIANT CONTEXT FILE====================
    def getVariantContexts(self):
        return [varcon for varcon in self.variantContexts.values()]

    # Returns a specified variant context
    def getVariantContext(self, contextId):
        if(contextId in self.variantContexts):
            return self.variantContexts[contextId]
        return None

    # Returns the acceptor context of the specified variant context
    def getAcceptorContext(self, contextId):
        if(contextId in self.variantContexts):
            return self.variantContexts[contextId].getAcceptorContext()
        return None

    # Returns the donor context of the specified variant context
    def getDonorContext(self, contextId):
        if(contextId in self.variantContexts):
            return self.variantContexts[contextId].getDonorContext()
        return None

    # Returns all variant context acceptor reads
    def getAllVariantContextAcceptorReads(self):
        acceptorReads = []
        for varcon in self.variantContexts.values():
            acceptorReads.extend(varcon.getVariantContextAcceptorReads())
        return acceptorReads

    # Returns all variant context donor reads
    def getAllVariantContextDonorReads(self):
        donorReads = []
        for varcon in self.variantContexts.values():
            donorReads.extend(varcon.getVariantContextDonorReads())
        return donorReads

    # Returns all variant context acceptor read ids
    def getAllVariantContextAcceptorReadIds(self):
        acceptorReadIds = []
        for varcon in self.variantContexts.values():
            acceptorReadIds.extend(varcon.getVariantContextAcceptorReadIds())
        return acceptorReadIds

    # Returns the variant context donor read ids
    def getAllVariantContextDonorReadIds(self):
        donorReadIds = []
        for varcon in self.variantContexts.values():
            donorReadIds.extend(varcon.getVariantContextDonorReadIds())
        return donorReadIds



    # ====================BASIC VARIANTCONTEXTFILE METHODS====================
    # Reads a provided variant context file and saves data according to set filters
    def readVariantContextFile(self, fileLoc, sampleFilter=None, varconFilter=None, chromFilter=None):
        try:
            with open(fileLoc, 'r') as vcFile:
                next(vcFile)    # Skip the header line
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


    # Returns whether something is in the filter or not
    def passesFilter(self, valToCheck, filterList):
        if(filterList is not None):
            return valToCheck in filterList
        return True



    # ====================VARIANT CONTEXT SET OPERATIONS (UNION, INTERSECT, DIFFERENCE)====================
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



    #====================METHODS TO OBTAIN VARIANT CONTEXT DATA BASED ON FILTERS====================
    # Returns a list/hashmap of VariantContextObjects
    def getVariantContexts(self, asList=False, varconFilter=None, sampleFilter=None, chromFilter=None):
        if(asList):
            return [x for x in self.variantContexts.values() if(self.passesFilter(x.getVariantContextId(), varconFilter) and self.passesFilter(x.getVariantContextSample(), sampleFilter) and self.passesFilter(x.getVariantContextChrom(), chromFilter))]
        return {k:v for k,v in self.variantContexts if(self.passesFilter(k, varconFilter) and self.passesFilter(v.getVariantContextSample(), sampleFilter) and self.passesFilter(v.getVariantContextChrom(), chromFilter))}



    # ====================METHODS TO ASSESS WHETHER A VARIANT IS IN AN EXISTING CONTEXT====================
    # Main method that returns whether a variant (SNP or indel).
    def variantIsInContext(self, variantType, searchChrom, searchStart, searchStop):
        if(variantType=='snp'):
            return self.snpVariantIsInContext(searchChrom, searchStart)
        if(variantType=='indel'):
            return self.indelVariantIsInContext(searchChrom, searchStart, searchStop)


    # Determines whether an SNP variant is located in an already existing variant context.
    def snpVariantIsInContext(self, varchrom, vcfVarPos):
        for varcon in self.variantContexts.values():
            if(varchrom == varcon.getVariantContextChrom()):
                if(vcfVarPos >= varcon.getVariantContextStart() and vcfVarPos <= varcon.getVariantContextEnd()):
                    return True
        return False


    # Determines whether an indel variant is located within an existing variant context (indelLeftPos and indelRightPos can be the used search window)
    def indelVariantIsInContext(self, indelChrom, indelLeftPos, indelRightPos):
        for varcon in self.variantContexts.values():
            if(indelChrom == varcon.getVariantContextChrom()):
                if(indelLeftPos <= varcon.getContextStart() and indelRightPos >= varcon.getContextStart()):
                    return True
                if(indelLeftPos <= varcon.getContextEnd() and indelRightPos >= varcon.getContextEnd()):
                    return True
                if(indelLeftPos >= varcon.getContextStart() and indelRightPos <= varcon.getContextEnd()):
                    return True
        return False



    # ====================METHODS TO ADD DATA/VARIANT CONTEXTS TO THE VARIANT CONTEXT FILE====================
    # Adds a variant context object
    def addVariantContext(self, varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconAReads, otherVarconDReads, varconALength=None, varconDLength=None):
        varconObj = VariantContext2(varconId, varconSample, varconChrom, varconOrigin, varconStart, varconEnd, varconAReads, otherVarconDReads, varconALength, varconDLength)
        self.variantContexts[varconId] = varconObj

    # Adds an acceptor context object to a variant context
    def setAcceptorContext(self, varconId, accContext):
        if(varconId in self.variantContexts):
            self.variantContexts[varconId].setAcceptorContext(accContext)

    # Add a newly created acceptor context to an existing variant context
    def addAcceptorContext(self, contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, acceptorReads):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].addAcceptorContext(contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, acceptorReads)

    # Sets the donor context object of a variant context
    def setDonorContext(self, varconId, donContext):
        if(varconId in self.variantContexts):
            self.variantContexts[varconId].setDonorContext(donContext)

    # Add a newly created donor context to an existing variant context
    def addDonorContext(self, contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, donorReads):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].addDonorContext(contextId, sampleId, contextChrom, contextOrigin, contextStart, contextEnd, donorReads)



    # ====================METHODS TO ADD UNMAPPED MATE IDS TO ACCEPTOR, DONOR AND VARIANT CONTEXT====================
    # Sets the unmapped mate read id for a specified acceptor context
    def setAcceptorContextUnmappedMateIds(self, contextId, mateIds):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].setAcceptorContextUnmappedMates(mateIds)

    # Sets the unmapped mate read id for a specified donor context
    def setDonorContextUnmappedMateIds(self, contextId, mateIds):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].setDonorContextUnmappedMates(mateIds)

    # Sets the acceptor unmapped mate ids for a specified variant context
    def setUnmappedAcceptorMateIds(self, contextId, mateIds):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].setUnmappedAcceptorMateIds(mateIds)

    # Sets the donor unmapped mate ids for a specified variant context
    def setUnmappedDonorMateIds(self, contextId, mateIds):
        if(contextId in self.variantContexts):
            self.variantContexts[contextId].setUnmappedDonorMateIds(mateIds)



    # ====================METHODS TO OBTAIN SOME STATISTICS ABOUT ALL THE CONTEXTS====================
    # Returns the average variant context length within this variant context file
    def getAverageVariantContextLength(self):
        return statistics.mean([varcon.getVariantContextLength() for varcon in self.variantContexts.values()])

    # Returns the median variant context length within this variant context file
    def getMedianVariantContextLength(self):
        return statistics.median([varcon.getVariantContextLength() for varcon in self.variantContexts.values()])

    # Returns the average number of variant context reads for this variant context file
    def getAverageVariantContextReads(self):
        return statistics.mean()

    # Returns the median number of variant context reads for this variant context file
    def getMedianVariantContextReads(self):
        return statistics.median()



    # ====================METHODS TO WRITE VARIANT CONTEXT DATA TO A FILE====================
    # Writes the variant context data to an output file
    def writeVariantContextFile(self, outFileLoc, sampleFilter=None, varconFilter=None, chromFilter=None):
        try:
            with open(outFileLoc, 'w') as varconOutFile:
                varconOutFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tAcceptorContextLength\tDonorContextLength\tAcceptorReads\tDonorReads\tADratio\tAcceptorReadsIds\tDonorReadIds\n")
                for varcon in self.variantContexts.values():
                    samplePass = self.passesFilter(varcon.getVariantContextSample(), sampleFilter)
                    varconPass = self.passesFilter(varcon.getVariantContextId(), varconFilter)
                    chromPass = self.passesFilter(varcon.getVariantContextChrom(), chromFilter)
                    if(samplePass and varconPass and chromPass):
                        varconOutFile.write(varcon.toString()+"\n")
        except IOError as ioe:
            self.vaseUtilLogger.warning("Could not write variant contexts to " +str(ioe.filename))

    # Writes the donor contexts used to construct the variant contexts
    def writeAcceptorContexFile(self, outFileLoc, sampleFilter=None, contextFilter=None, chromFilter=None):
        try:
            with open(outFileLoc, 'w') as varconOutFile:
                varconOutFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tNumOfReads\tReadIds\n")
                for varcon in self.variantContexts.values():
                    samplePass = self.passesFilter(varcon.getVariantContextSample(), sampleFilter)
                    varconPass = self.passesFilter(varcon.getVariantContextId(), varconFilter)
                    chromPass = self.passesFilter(varcon.getVariantContextChrom(), chromFilter)
                    if(samplePass and varconPass and chromPass):
                        varconOutFile.write(varcon.getAcceptorContext().toString()+"\n")
        except IOError as ioe:
            self.vaseUtilLogger.warning("Could not write acceptor contexts to " +str(ioe.filename))

    # Writes the acceptor cotnexts used to construct the variant contexts
    def writeDonorContexFile(self, outFileLoc, sampleFilter=None, contextFilter=None, chromFilter=None):
        try:
            with open(outFileLoc, 'w') as varconOutFile:
                varconOutFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tNumOfReads\tReadIds\n")
                for varcon in self.variantContexts.values():
                    samplePass = self.passesFilter(varcon.getVariantContextSample(), sampleFilter)
                    varconPass = self.passesFilter(varcon.getVariantContextId(), varconFilter)
                    chromPass = self.passesFilter(varcon.getVariantContextChrom(), chromFilter)
                    if(samplePass and varconPass and chromPass):
                        varconOutFile.write(varcon.getDonorContext().toString()+"\n")
        except IOError as ioe:
            self.vaseUtilLogger.warning("Could not write donor contexts to " +str(ioe.filename))


    # Writes some statistics about the acceptor and donor reads identified for each variant context
    def writeVariantContextStats(self, statsOutLoc):
        try:
            with open(statsOutLoc, 'w') as varconStatsFile:
                for varcon in self.variantContexts.values():
                    varconStatsFile.write(varcon.toStatisticsString()+"\n")
        except IOError as ioe:
            self.vaseLogger.critical("Coud not write variant context statistics to " +str(statsOutLoc))

    # Writes some statistics about the acceptor and donor reads identified for each variant context
    def writeAcceptorContextStats(self, statsOutLoc):
        try:
            with open(statsOutLoc, 'w') as varconStatsFile:
                for varcon in self.variantContexts.values():
                    varconStatsFile.write(varcon.getAcceptorContext().toStatisticsString()+"\n")
        except IOError as ioe:
            self.vaseLogger.critical("Coud not write acceptor context statistics to " +str(statsOutLoc))

    # Writes some statistics about the acceptor and donor reads identified for each variant context
    def writeDonorContextStats(self, statsOutLoc):
        try:
            with open(statsOutLoc, 'w') as varconStatsFile:
                for varcon in self.variantContexts.values():
                    varconStatsFile.write(varcon.getDonorContext().toStatisticsString()+"\n")
        except IOError as ioe:
            self.vaseLogger.critical("Coud not write donor context statistics to " +str(statsOutLoc))


    # Writes the left and right positions to the output file. Left pos for R1 and right pos for R2
    def writeLeftRightPositions(self, typeToWrite, outFileLoc):
        try:
            with open(outFileLoc, 'w') as lrpof:
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variantContexts.values():
                    leftPositions, rightPositions = [], []
                    if(typeToWrite=='acceptor'):
                        leftPositions = varcon.getAcceptorReadRefStarts()
                        rightPositions = varcon.getAcceptorReadRefEnds()
                    if(typeToWrite=='donor'):
                        leftPositions = varcon.getDonorReadRefStarts()
                        rightPositions = varcon.getDonorReadRefEnds()
                    lrpof.write(str(contextId) +"\t"+ ','.join(leftPositions) +"\t"+ ','.join(rightPositions))
        except IOError as ioe:
            self.vaseLogger.warning("Could not write read left positions to output file " +str(outFileLoc))

    # Writes the acceptor context left and right positions to the output file. Left pos for R1 and right pos for R2
    def writeAcceptorLeftRightPositions(self, outFileLoc):
        try:
            with open(outFileLoc, 'w') as lrpof:
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variantContexts.values():
                    leftPositions = varcon.getAcceptorContext().getContextBamReadLeftPositions()
                    rightPositions = varcon.getAcceptorContext().getContextBamReadRightPositions()
                    lrpof.write(str(contextId) +"\t"+ ','.join(leftPositions) +"\t"+ ','.join(rightPositions))
        except IOError as ioe:
            self.vaseLogger.warning("Could not write read left positions to output file " +str(outFileLoc))

    # Writes the left and right positions to the output file. Left pos for R1 and right pos for R2
    def writeDonorLeftRightPositions(self, outFileLoc):
        try:
            with open(outFileLoc, 'w') as lrpof:
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variantContexts.values():
                    leftPositions = varcon.getAcceptorContext().getContextBamReadLeftPositions()
                    rightPositions = varcon.getAcceptorContext().getContextBamReadRightPositions()
                    lrpof.write(str(contextId) +"\t"+ ','.join(leftPositions) +"\t"+ ','.join(rightPositions))
        except IOError as ioe:
            self.vaseLogger.warning("Could not write read left positions to output file " +str(outFileLoc))


    # Writes the identifiers of reads that have unmapped mates per sample to a file. Samples are all donors and the ?template?.
    def writeReadsWithUnmappedMate(self, typeToWrite, umFileLoc):
        try:
            with open(umFileLoc, 'w') as umFile:
                umFile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variantContexts.values():
                    if(typeToWrite=='acceptor'):
                        umFile.write(varcon.getVariantContextId() +"\t"+ str(varcon.getVariantContextSample()) +"\t"+ ';'.join(varcon.getUnmappedAcceptorMateIds()))
                    if(typeToWrite=='donor'):
                        umFile.write(varcon.getVariantContextId() +"\t"+ str(varcon.getVariantContextSample()) +"\t"+ ';'.join(varcon.getUnmappedDonorMateIds()))
        except IOError:
            self.vaseLogger.warning("Could not write read identifiers of reads with unmapped mates to " +str(umFileLoc))

    # Writes the unmapped mate id of the acceptor context
    def writeAcceptorUnmappedMates(self, umFileLoc):
        try:
            with open(umFileLoc, 'w') as umFile:
                umFile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variantContexts.values():
                    acccon = varcon.getAcceptorContext()
                    umFile.write(str(acccon.getContextId())+ "\t" +str(acccon.getSampleId())+ "\t" +';'.join(acccon.getUnmappedReadMateIds()))
        except IOError:
            self.vaseLogger.warning("Could not write read identifiers of reads with unmapped mates to " +str(umFileLoc))

    # Writes the unmapped mate id of the acceptor context
    def writeDonorUnmappedMates(self, outFileLoc):
        try:
            with open() as umFile:
                umFile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variantContexts.values():
                    doncon = varcon.getDonorContext()
                    umFile.write(str(doncon.getContextId())+ "\t" +str(doncon.getSampleId())+ "\t" +';'.join(doncon.getUnmappedReadMateIds()))
        except IOError:
            self.vaseLogger.warning("Could not write read identifiers of reads with unmapped mates to " +str(umFileLoc))
