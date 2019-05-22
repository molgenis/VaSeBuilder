#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import statistics
import pysam

# Import VaSe specific classes.
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile


class VaSeBuilder:
    # Constructor that saves the identifier, date and time of the current
    def __init__(self, vaseId):
        self.vaseLogger = logging.getLogger("VaSe_Logger")
        self.creationId = str(vaseId)
        self.creationDate = datetime.now().date()
        self.creationTime = datetime.now().time()

        # Create the bookkeeping variables used for saving the variant
        # contexts.
        # VariantContextFile that saves the acceptor, donor and variant
        # contexts and their associated data.
        self.contexts = VariantContextFile()
        self.vaseLogger.info(f"VaSeBuilder: "
                             + str(self.creationId) + " ; "
                             + str(self.creationDate) + " ; "
                             + str(self.creationTime))

    # ===METHODS TO PRODUCE VARIANT CONTEXTS===================================
    # Creates the new FastQ validation dataset by replacing NIST reads
    # containing a VCF variant with BAM reads from patients.  Returns
    # true at the end to indicate the process is done.
    def buildValidationSet(self, vcfBamLinkMap,
                           vcfSampleMap, bamSampleMap,
                           acceptorBamLoc,
                           fastqFPath, fastqRPath,
                           outPath,
                           fastqOutPath,
                           varConOutPath):
        self.vaseLogger.info("Start building the validation set")
        donorVcfsUsed, donorBamsUsed = [], []

        try:
            acceptorBamFile = pysam.AlignmentFile(acceptorBamLoc, "rb")

            # Iterate over the samples to use for building the validation set.
            for sampleId in vcfSampleMap:
                self.vaseLogger.debug(f"Processing data for sample {sampleId}")

                # Check in the VCF BAM link map that there is a BAM file
                # for the sample as well.
                try:
                    vcfFile = pysam.VariantFile(vcfSampleMap[sampleId], "r")
                    self.vaseLogger.debug("Opened VCF file "
                                          f"{vcfSampleMap[sampleId]}")
                    bamFile = pysam.AlignmentFile(bamSampleMap[sampleId], "rb")
                    self.vaseLogger.debug("Opened BAM file "
                                          f"{bamSampleMap[sampleId]}")

                    # Loop over the variants in the VCF file. Prior to
                    # identifying the BAM reads, it is first checked
                    # whether the variant is in a previously established
                    for vcfVar in vcfFile.fetch():
                        variantId = self.getVcfVariantId(vcfVar)
                        acceptorUnmapped = []
                        donorUnmapped = []
                        varconUnmapped_A = []
                        varconUnmapped_D = []
                        self.vaseLogger.debug("Searching BAM reads for "
                                              f"variant {variantId}")

                        # Determine the search window before gathering
                        # reads overlapping with the VCF variant.
                        variantType = self.determineVariantType(vcfVar.ref,
                                                                vcfVar.alts)
                        searchWindow = self.determineReadSearchWindow(
                                variantType,
                                vcfVar
                                )

                        # Get the BAM reads for the variant and
                        # determine the variant context.
                        if (not self.contexts.variantIsInContext(
                                variantType,
                                vcfVar.chrom,
                                searchWindow[0],
                                searchWindow[1]
                                )):
                            try:
                                # Gather acceptor reads and their mates
                                # overlapping with the variant and determine
                                # the acceptor context.
                                self.vaseLogger.debug(
                                        "Determine acceptor BAM reads for "
                                        f"variant {variantId}"
                                        )
                                # Obtain all acceptor BAM reads containing the
                                # VCF variant and their read mate.
                                acceptorContextReads = self.getVariantReads(
                                        variantId, vcfVar.chrom,
                                        searchWindow[0], searchWindow[1],
                                        acceptorBamFile, True,
                                        acceptorUnmapped
                                        )
                                self.vaseLogger.debug(
                                        "Determine acceptor context for "
                                        f"variant {variantId}"
                                        )
                                # Determine the acceptor variant context based
                                # on the reads overlapping the variant.
                                acceptorContext = self.determineContext(
                                        acceptorContextReads,
                                        vcfVar.pos
                                        )

                                # Gather donor reads and their mates
                                # overlapping with the variant and determine
                                # the donor context.
                                self.vaseLogger.debug(
                                        "Search donor BAM reads for variant "
                                        f"{variantId}"
                                        )
                                # Obtain all donot BAM reads containing the VCF
                                # variant and their read mate.
                                donorContextReads = self.getVariantReads(
                                        variantId, vcfVar.chrom,
                                        searchWindow[0], searchWindow[1],
                                        bamFile, True,
                                        donorUnmapped
                                        )
                                self.vaseLogger.debug(
                                        "Determine donor context for variant "
                                        f"{variantId}"
                                        )
                                # Determine the donor variant context based on
                                # the reads overlapping the variant.
                                donorContext = self.determineContext(
                                        donorContextReads,
                                        vcfVar.pos
                                        )

                                # Determine the ultimate variant context and
                                # obtain the overlapping acceptor and donor
                                # reads.
                                variantContext = self.determineLargestContext(
                                        vcfVar.pos,
                                        acceptorContext,
                                        donorContext
                                        )
                                # Obtain all acceptor reads overlapping with
                                # the combined variant context and their mates.
                                variantContextAcceptorReads = self.getVariantReads(
                                        variantId,
                                        variantContext[0],
                                        variantContext[2],
                                        variantContext[3],
                                        acceptorBamFile, True,
                                        varconUnmapped_A
                                        )
                                # Obtain all donor reads overlapping with the
                                # combined variant context and their mates.
                                variantContextDonorReads = self.getVariantReads(
                                        variantId,
                                        variantContext[0],
                                        variantContext[2],
                                        variantContext[3],
                                        bamFile, True,
                                        varconUnmapped_D
                                        )

                                # Check whether reads were found in both
                                # acceptor and donor.  Only then save the
                                # results.
                                if ((len(donorContextReads) > 0)
                                   and (len(acceptorContextReads) > 0)):
                                    self.contexts.addVariantContext(
                                            variantId,
                                            sampleId,
                                            variantContext[0],
                                            variantContext[1],
                                            variantContext[2],
                                            variantContext[3],
                                            variantContextAcceptorReads,
                                            variantContextDonorReads
                                            )
                                    self.contexts.addAcceptorContext(
                                            variantId,
                                            sampleId,
                                            acceptorContext[0],
                                            acceptorContext[1],
                                            acceptorContext[2],
                                            acceptorContext[3],
                                            acceptorContextReads
                                            )
                                    self.contexts.addDonorContext(
                                            variantId,
                                            sampleId,
                                            donorContext[0],
                                            donorContext[1],
                                            donorContext[2],
                                            donorContext[3],
                                            donorContextReads
                                            )

                                    # Add the read identifiers of reads with
                                    # an unmapped mate.
                                    self.contexts.setUnmappedAcceptorMateIds(
                                            variantId,
                                            varconUnmapped_A
                                            )
                                    self.contexts.setUnmappedDonorMateIds(
                                            variantId,
                                            varconUnmapped_D
                                            )
                                    self.contexts.setAcceptorContextUnmappedMateIds(
                                            variantId,
                                            acceptorUnmapped
                                            )
                                    self.contexts.setDonorContextUnmappedMateIds(
                                            variantId,
                                            donorUnmapped
                                            )
                                else:
                                    self.vaseLogger.debug(
                                            "No donor and/or acceptor BAM "
                                            "reads found for variant "
                                            f"{variantId}"
                                            )

                            except IOError as ioe:
                                self.vaseLogger.warning(
                                        "Could not obtain BAM reads from "
                                        f"{bamSampleMap[sampleId]}"
                                        )
                        else:
                            self.vaseLogger.debug(
                                    f"VCF variant {variantId} is located in "
                                    "an already existing variant context"
                                    )
                    bamFile.close()
                    vcfFile.close()

                    donorVcfsUsed.append(vcfSampleMap[sampleId])
                    donorBamsUsed.append(bamSampleMap[sampleId])
                except IOError as ioe:
                    self.vaseLogger.warning("Could not establish data for "
                                            f"{sampleId}")

            # Write the variant context data and used donor VCFs/BAMs
            # to output files.
            self.vaseLogger.info(
                    f"Writing combined variant contexts to {varConOutPath}"
                    )
            # Writes the variant contexts to file.
            self.contexts.writeVariantContextFile(varConOutPath)

            self.vaseLogger.info(
                    "Writing combined variant context statistics to "
                    f"{outPath}/varconstats.txt"
                    )
            # Writes the variant context statistics file.
            self.contexts.writeVariantContextStats(
                    f"{outPath}/varconstats.txt"
                    )

            self.vaseLogger.info(
                    "Write the used donor VCF files per sample to "
                    f"{outPath}/donorvcfs.txt"
                    )
            self.writeUsedDonorFiles(
                    f"{outPath}/donorvcfs.txt",
                    vcfSampleMap,
                    donorVcfsUsed
                    )

            self.vaseLogger.info(
                    "Write the used donor BAM files per sample to "
                    f"{outPath}/donorbams.txt"
                    )
            self.writeUsedDonorFiles(
                    f"{outPath}/donorbams.txt",
                    bamSampleMap,
                    donorBamsUsed
                    )

            # Checks whether the program is running on debug.  If so,
            # write some extra output files.
            if (self.vaseLogger.getEffectiveLevel() == 10):
                self.writeOptionalOutputFiles(outPath, self.contexts)

            # Obtain a list of acceptor reads to skip when iterating
            # over the acceptor FastQ.
            # Set up a list of all acceptor reads to skip.
            acceptorReadsToSkip = list(set(
                    self.contexts.getAllVariantContextAcceptorReadIds()
                    ))
            # Sets up a list.
            donorReads = self.contexts.getAllVariantContextDonorReads()

            # Make the new FastQ files that can be used to run in the
            # NGS_DNA pipeline along real sample data.
            self.vaseLogger.info("Start writing the R1 FastQ files")
            # Build the R1 fastq file.
            self.buildFastQ(fastqFPath, acceptorReadsToSkip,
                            donorReads, "F", fastqOutPath)
            self.vaseLogger.info("Wrote all R1 FastQ files")

            self.vaseLogger.info("Start writing the R2 FastQ files")
            # Build the R2 fastq file.
            self.buildFastQ(fastqRPath, acceptorReadsToSkip,
                            donorReads, "R", fastqOutPath)
            self.vaseLogger.info("Wrote all R2 FastQ files")

            self.vaseLogger.info("Finished building the validation set")
            acceptorBamFile.close()

        except IOError as ioe:
            self.vaseLogger.critical("Could not open acceptor BAM file")
            exit()

    # Returns whether a variant is a SNP or indel.
    def determineVariantType(self, vcfVariantRef, vcfVariantAlts):
        # Determine the maximum reference allele length.
        maxRefLength = max([len(x) for x in vcfVariantRef.split(",")])
        # Determine the maximum alternative allele length.
        maxAltLength = max([len(x) for x in vcfVariantAlts])

        # Check based on the reference and alternative lengths whether
        # the variant is a SNP or indel.
        if (maxRefLength == 1 and maxAltLength == 1):
            return "snp"
        elif (maxRefLength > 1 or maxAltLength > 1):
            return "indel"
        return "?"

    # Determines the start and end positions to use for searching reads
    # overlapping with the variant.
    def determineReadSearchWindow(self, variantType, vcfVariant):
        if (variantType == "snp"):
            return [vcfVariant.pos - 1, vcfVariant.pos + 1]
        elif (variantType == "indel"):
            return self.determineIndelReadRange(vcfVariant.pos,
                                                vcfVariant.ref,
                                                vcfVariant.alts)
        return [-1, -1]

    # Returns the search start and stop to use for searching BAM reads
    # overlapping with the range of the indel.
    def determineIndelReadRange(self, variantPos, variantRef, variantAlts):
        searchStart = variantPos
        searchStop = variantPos + max(
                max([len(x) for x in variantRef.split(",")]),
                max([len(x) for x in variantAlts])
                )
        return [searchStart, searchStop]

    # Returns the BAM reads containing the specific vcf variant as well
    # as their read mate.
    def getVariantReads(self, contextid, vcfVariantChr,
                        varStartPos, varEndPos,
                        bamFile,
                        writeUnm=False, uMateList=None):
        # Obtain all the variant reads overlapping with the variant and
        # their mate reads.
        variantReads = []

        for vread in bamFile.fetch(vcfVariantChr, varStartPos, varEndPos):
            variantReads.append(DonorBamRead(
                    vread.query_name,
                    self.getReadPairNum(vread),
                    vread.reference_name,
                    vread.reference_start,
                    vread.infer_read_length(),
                    vread.get_forward_sequence(),
                    "".join([chr((x + 33))
                             for x in vread.get_forward_qualities()]),
                    vread.mapping_quality
                    ))

            # Try to obtain the reads mate as well.
            try:
                vmread = bamFile.mate(vread)
                variantReads.append(DonorBamRead(
                        vmread.query_name,
                        self.getReadPairNum(vmread),
                        vmread.reference_name,
                        vmread.reference_start,
                        vmread.infer_read_length(),
                        vmread.get_forward_sequence(),
                        "".join([chr((x+33))
                                 for x in vmread.get_forward_qualities()]),
                        vmread.mapping_quality
                        ))

            except ValueError as pve:
                self.vaseLogger.debug("Could not find mate for "
                                      f"{vread.query_name} ; "
                                      "mate is likely unmapped.")
                if (writeUnm):
                    uMateList.append(vread.query_name)

        # Make sure the list only contains each BAM read once (if a read
        # and mate both overlap with a variant, they have been added
        # twice to the list).
        # OLD: variantReads = list(set(variantReads))
        uniq_variantReads = []
        for fetched in variantReads:
            if ((fetched.getBamReadId(), fetched.getBamReadPairNumber())
               not in [(y.getBamReadId(), y.getBamReadPairNumber())
                       for y in uniq_variantReads]):
                uniq_variantReads.append(fetched)
        variantReads = uniq_variantReads

        # Filter to keep only read pairs.
        variantReads = self.filterVariantReads(variantReads)

        self.vaseLogger.debug(f"Found a total of {len(variantReads)} "
                              "BAM reads.")
        return variantReads

    # Filters the donor reads to keep only reads that occur twice.
    def filterVariantReads(self, bamReads):
        return [bread for bread in bamReads
                if (self.readOccurence(bread.getBamReadId(), bamReads) == 2)]

    # Returns the number of occurences of a certain read in the list of
    # BAM reads (should be two ideally).
    def readOccurence(self, readId, readList):
        return sum([bamRead.getBamReadId() == readId for bamRead in readList])

    # Determines the start and stops of the variant context (please see
    # the documentation for more information).
    def determineContext(self, contextReads, contextOrigin):
        # Check whether there are reads to determine the context for.
        if (len(contextReads) > 0):
            contextChrom = contextReads[0].getBamReadChrom()
            contextStart = min([conread.getBamReadRefPos()
                                for conread in contextReads])
            contextEnd = max([conread.getBamReadRefEnd()
                              for conread in contextReads])
            self.vaseLogger.debug(f"Context is {contextChrom}, {contextStart},"
                                  f" {contextEnd}")
            return [contextChrom, contextOrigin, contextStart, contextEnd]
        return []

    # Determines the size of the variant context based on both the
    # acceptor and donor reads.
    def determineLargestContext(self, contextOrigin, acceptorContext,
                                donorContext):
        largestContext = [acceptorContext[0]]
        largestContext.append(contextOrigin)
        # Determine and save the smallest context start.
        largestContext.append(min(acceptorContext[2], donorContext[2]))
        # Determine and save the largest context end.
        largestContext.append(max(acceptorContext[3], donorContext[3]))
        return largestContext

    # ===METHODS TO PRODUCE THE VALIDATION FASTQ FILES=========================
    # Will build the R1/R2 VaSe fastq files.
    def buildFastQ(self, acceptorFqFilePaths, acceptorReadsToSkip,
                   donorContextReadMap, forwardOrReverse, vaseFqOutPath):
        writeDonor = False

        # Iterate over the R1/R2 fastq in files to use as templates for the
        for x in range(0, len(acceptorFqFilePaths)):
            if (x == len(acceptorFqFilePaths)-1):
                writeDonor = True
                self.vaseLogger.debug("Donor reads will be added the current "
                                      "VaSe fastQ out file.")

            # Write the new VaSe FastQ file.
            vaseFqOutName = self.setFastqOutPath(vaseFqOutPath,
                                                 forwardOrReverse,
                                                 x + 1)
            self.vaseLogger.debug(f"Set FastQ output path to: {vaseFqOutName}")
            self.writeVaSeFastQ(acceptorFqFilePaths[x], vaseFqOutName,
                                acceptorReadsToSkip, donorContextReadMap,
                                forwardOrReverse, writeDonor)

    # Builds a new FastQ file to be used for validation.
    def writeVaSeFastQ(self, acceptorFastqIn, fastqOutPath,
                       acceptorReadsToSkip, donorBamReadData,
                       fR, writeDonorData=False):
        try:
            fqFile = io.BufferedWriter(gzip.open(fastqOutPath, "wb"))
            self.vaseLogger.debug(f"Opened template FastQ: {acceptorFastqIn}")

            # Open the template fastq and write filtered data to a new
            # fastq.gz file.
            gzFile = io.BufferedReader(gzip.open(acceptorFastqIn, "rb"))
            for fileLine in gzFile:

                # Check if we are located at a read identifier.
                if (fileLine.startswith(b"@")):
                    if (fileLine.decode("utf-8").strip()[1:] not in acceptorReadsToSkip):
                        fqFile.write(fileLine)
                        fqFile.write(next(gzFile))
                        fqFile.write(next(gzFile))
                        fqFile.write(next(gzFile))
            gzFile.close()

            # Add the patient BAM reads containing a VCF variant to the
            # new FastQ file.
            if (writeDonorData):
                donorBamReadData.sort(key=lambda x: x.getBamReadId(),
                                      reverse=False)
                for bamRead in donorBamReadData:
                    # Check if the BAM read is R1 or R2.
                    if (self.isRequiredRead(bamRead, fR)):
                        fqFile.write(bamRead.getAsFastQSeq().encode("utf-8"))
            fqFile.flush()
            fqFile.close()

        except IOError as ioe:
            if (ioe.filename == acceptorFastqIn):
                self.vaseLogger.critical("The supplied template FastQ file "
                                         "could not be found.")
            if (ioe.filename == fastqOutPath):
                self.vaseLogger.critical("A FastQ file could not be written "
                                         "to the provided output location.")
            exit()

    # Checks if a read is read 1 (R1) or read 2 (R2).
    def isRequiredRead(self, bamRead, fR):
        if (fR == "F"):
            return bamRead.isRead1()
        return bamRead.isRead2()

    # Returns the name for the fastq out file.
    def setFastqOutPath(self, outPath, fR, lNum):
        if (fR == "F"):
            return (f"{outPath}_{datetime.now().date()}_L{lNum}_R1.fastq.gz")
        return (f"{outPath}_{datetime.now().date()}_L{lNum}_R2.fastq.gz")

    # ===METHODS TO OBTAIN SOME DATA OF THE VASEBUILDER OBJECT=================
    # Returns the identifier of the current VaSeBuilder object.
    def getCreationId(self):
        return self.creationId

    # Returns the date the current VaSeBuilder object has been made.
    def getCreationDate(self):
        return self.creationDate

    # Returns the time the current VaSeBuilder object has been made.
    def getCreationTime(self):
        return self.creationTime

    # Returns all variant contexts.
    def getVariantContexts(self):
        return self.variantContextMap

    # Returns a specified acceptor context.
    def getAcceptorContext(self, contextId):
        return self.acceptorContexts.getVariantContext(contextId)

    # Returns a specified donor context.
    def getDonorContext(self, contextId):
        return self.donorContexts.getVariantContext(contextId)

    # Returns the context start and stop for a specified VCF variant.
    def getVariantContext(self, contextId):
        return self.variantContexts.getVariantContext(contextId)

    # Returns an identifier for a VCF variant.
    # If the identifier is "." then one will be constructed as
    # "chrom_pos".
    def getVcfVariantId(self, vcfVariant):
        return (f"{vcfVariant.chrom}_{vcfVariant.pos}")

    # Returns whether the read is the first or second read in a pair.
    def getReadPairNum(self, pysamBamRead):
        if (pysamBamRead.is_read1):
            return "1"
        return "2"

    # ===METHODS TO WRITE OUTPUT FILES=========================================
    # Writes the used donor vcf files to a file
    def writeUsedDonorFiles(self, outLocFile, fileSampleMap,
                            listOfUsedDonorFiles):
        try:
            with open(outLocFile, "w") as outFile:
                outFile.write("#SampleId\tDonorFile\n")
                for sampleid, sampleFile in fileSampleMap.items():
                    if (sampleFile in listOfUsedDonorFiles):
                        outFile.write(f"{sampleid}\t{sampleFile}\n")
        except IOError as ioe:
            self.vaseLogger.critical("Could not write used donor files to "
                                     f"{outLocFile}")

    # Writes the optional output files (when logger is set to DEBUG log
    # level).
    def writeOptionalOutputFiles(self, outPath, contextFile):
        # Write the optional acceptor context files; acceptor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaseLogger.debug("Writing acceptor contexts to "
                              f"{outPath}/acceptorcontexts.txt")
        contextFile.writeAcceptorContextFile(f"{outPath}/acceptorcontexts.txt")

        self.vaseLogger.debug("Writing acceptor context statistics to "
                              f"{outPath}/acceptorcontextstats.txt")
        contextFile.writeAcceptorContextStats(
                              f"{outPath}/acceptorcontextstats.txt"
                              )

        self.vaseLogger.debug("Writing acceptor context read identifiers with "
                              "unmapped mates to"
                              f"mates to {outPath}/acceptor_unmapped.txt")
        contextFile.writeAcceptorUnmappedMates(
                              f"{outPath}/acceptor_unmapped.txt"
                              )

        self.vaseLogger.debug("Writing left and right most read positions of "
                              "each acceptor context to "
                              f"{outPath}/acceptor_positions.txt")
        contextFile.writeAcceptorLeftRightPositions(
                              f"{outPath}/acceptor_positions.txt"
                              )

        # Write the optional donor context files; donor contexts,
        # read ids with unmapped mate and left/right positions.
        self.vaseLogger.debug("Writing donor contexts to "
                              f"{outPath}/donorcontexs.txt")
        contextFile.writeDonorContextFile(f"{outPath}/donorcontexts.txt")

        self.vaseLogger.debug("Writing donor context statistics to "
                              f"{outPath}/donorcontextstats.txt")
        contextFile.writeDonorContextStats(f"{outPath}/donorcontextstats.txt")

        self.vaseLogger.debug("Writing donor context read identifiers "
                              "with unmapped mates to "
                              f"{outPath}/donor_unmapped.txt")
        contextFile.writeDonorUnmappedMates(f"{outPath}/donor_unmapped.txt")

        self.vaseLogger.debug("Writing left and right most read positions "
                              "of each donor context to "
                              f"{outPath}/donor_positions.txt")
        contextFile.writeDonorLeftRightPositions(
                              f"{outPath}/donor_positions.txt"
                              )

        # Write the optional variant context files; acceptor & donor
        # unmapped mates and left/right positions.
        self.vaseLogger.debug("Writing variant context acceptor read "
                              "identifiers with unmapped mates to "
                              f"{outPath}/varcon_unmapped_acceptor.txt")
        contextFile.writeReadsWithUnmappedMate(
                              "acceptor",
                              f"{outPath}/varcon_unmapped_acceptor.txt"
                              )

        self.vaseLogger.debug("Writing variant context donor read identifiers "
                              "with unmapped mates to "
                              f"{outPath}/varcon_unmapped_donor.txt")
        contextFile.writeReadsWithUnmappedMate(
                              "donor",
                              f"{outPath}/varcon_unmapped_donor.txt"
                              )

        self.vaseLogger.debug("Writing variant context left and right most "
                              "read positions of acceptor reads to "
                              f"{outPath}/varcon_positions_acceptor.txt")
        contextFile.writeLeftRightPositions(
                              "acceptor",
                              f"{outPath}/varcon_positions_acceptor.txt"
                              )

        self.vaseLogger.debug("Writing variant context left and right most "
                              "read positions of donor reads to "
                              f"{outPath}/varcon_positions_donor.txt")
        contextFile.writeLeftRightPositions(
                              "donor",
                              f"{outPath}/varcon_positions_donor.txt"
                              )
