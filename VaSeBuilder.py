#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import statistics
import pysam
from DonorBamRead import DonorBamRead

class VaSeBuilder:
	# Constructor that saves the identifier, date and time of the current 
	def __init__(self, vaseId):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.creationId = str(vaseId)
		self.creationDate = datetime.now().date()
		self.creationTime = datetime.now().time()
		
		# Create the bookkeeping variables used for saving the variant contexts
		self.acceptorContextMap = {}	# Saves the acceptor variant context as {contextId => [chrom, origin, start, end]}
		self.donorContextMap = {}	# Saves the donor variant contexts as {contextId => [chrom, origin, start, end]}
		self.variantContextMap = {}	# Saves the largest variant contexts (the combination of acceptor and donor context) as {contextId => [chrom, origin, start, end]}
		
		# Create the bookkeeping variables for saving the read identifiers associated with variant contexts.
		self.acceptorContextReadMap = {}	# Saves the acceptor BAM reads per variant context as {variantId => [abamReadId1, abamReadId2, ..., abamReadIdn]}
		self.donorContextReadMap = {}	# Saves the donor BAM reads per variant context as {variantId => [dbamReadId1, dbamReadId2, ..., dbamReadIdn]}
		self.variantContextAcceptorReads = {}	# Saves the acceptor BAM read identifiers overlapping with a combined variant context as {contetId => [abamRead1, abamRead2, ..., abamReadn]}
		self.variantContextDonorReads = {}	# Saves the donor BAM read identifiers overlapping with a combined variant context as {contextId => [dbamRead1, dbamRead2, ..., dbamReadn]}
		
		# Create the bookkeeping variables for saving read identifiers with unmapped mates per variant context.
		self.variantContextUnmappedAcceptor = {}	# Saves the acceptor BAM read identifiers that have unmapped mates for the combined variant context as {contextId => [aubamRead1, auBamRead2, ..., auBamReadn]}
		self.variantContextUnmappedDonor = {}	# Saves the donor BAM read identifiers that have unmapped mates for the combined variant context as {contextId => [duBamRead1, duBamRead2, ..., duBamReadn]}
		self.acceptorUnmappedMateMap = {}	# Saves the read identifers of acceptor reads with unmapped mates per variant context as {contextId => [aubamRead1, aubamRead2, ..., aubamReadn]} 
		self.donorUnmappedMateMap = {}	# Saves the read identifiers of donor reads with unmapped mates per variant context as {contextId => [duBamRead1, duBamRead2, ..., duBamReadn]}
		
		self.vaseLogger.info("VaSeBuilder: " +str(self.creationId)+ " ; " +str(self.creationDate)+ " ; " +str(self.creationTime))
	
	
	# Creates the new FastQ validation dataset by replacing NIST reads containing a VCF variant with BAM reads from patients. Returns true at the end to indicate the process is done.
	def buildValidationSet(self, vcfBamLinkMap, vcfSampleMap, bamSampleMap, nistBamLoc, fastqFPath, fastqRPath, outPath, fastqOutPath, varConOutPath, varBreadOutPath, nistBreadOutPath):
		self.vaseLogger.info("Start building the validation set")
		variantSampleMap = {}	# Will link the variants to their samples. (Is used when writing the variant data to files)
		donorVcfsUsed, donorBamsUsed = [], []
		
		try:
			acceptorBamFile = pysam.AlignmentFile(nistBamLoc, 'rb')
			
			# Iterate over the samples to use for building the validation set.
			for sampleId in vcfSampleMap:
				self.vaseLogger.debug("Processing data for sample " +sampleId)
				
				# Check in the VCF BAM link map that there is a BAM file for the sample as well.
				try:
					vcfFile = pysam.VariantFile(vcfSampleMap[sampleId], 'r')
					self.vaseLogger.debug("Opened VCF file " +vcfSampleMap[sampleId]+ ".")
					
					bamFile = pysam.AlignmentFile(bamSampleMap[sampleId], 'rb')
					self.vaseLogger.debug("Opened BAM file " +bamSampleMap[sampleId])
					
					
					# Loop over the variants in the VCF file. Prior to identifying the BAM reads, it is first checked whether the variant is in a previously established 
					for vcfVar in vcfFile.fetch():
						variantId = self.getVcfVariantId(vcfVar)
						self.vaseLogger.debug("Searching BAM reads for variant " + str(variantId))
						
						# Get the BAM reads fr the variant and determine the variant context.
						if(not self.isInContext(vcfVar.chrom, vcfVar.pos)):
							try:
								# Determine the variant type and search window before gathering reads overlapping with the VCF variant
								variantType = self.determineVariantType(vcfVar.ref, vcfVar.alts)
								searchWindow = self.determineReadSearchWindow(variantType, vcfVar)
								
								
								self.vaseLogger.debug("Determine acceptor BAM reads for variant " +str(variantId))
								acceptorVariantReads = self.getVariantReads(variantId, vcfVar.chrom, searchWindow[0], searchWindow[1], acceptorBamFile, True, self.acceptorUnmappedMateMap)	# Obtain all NIST BAM reads containing the VCF variant and their read mate.
								
								self.vaseLogger.debug("Search donor BAM reads for variant " +str(variantId))
								donorVariantReads = self.getVariantReads(variantId, vcfVar.chrom, searchWindow[0], searchWindow[1], bamFile, True, self.donorUnmappedMateMap)	# Obtain all patient BAM reads containing the VCF variant and their read mate.
								
								# Determine the acceptor variant context
								self.vaseLogger.debug("Determine acceptor context for variant " +str(variantId))
								acceptorVariantContext = self.determineContext(acceptorVariantReads, vcfVar.pos)
								
								# Determine the donor variant context
								self.vaseLogger.debug("Determine donor context for variant " +str(variantId))
								donorVariantContext = self.determineContext(donorVariantReads, vcfVar.pos)	# Save the context start and stop for the variant in the VCF variant context map.
								
								# Determine the ultimate context and obtain acceptor and donor reads 
								variantContext = self.determineLargestContext(vcfVar.pos, acceptorVariantContext, donorVariantContext)
								
								
								# Determine the acceptor and donor BAM reads overlapping with the combined variant context
								acceptorContextReads = self.getVariantReads(variantId, variantContext[0], variantContext[2], variantContext[3], acceptorBamFile, True, self.variantContextUnmappedAcceptor)	# Obtain all acceptor reads overlapping with the combined variant context and their mates.
								donorContextReads = self.getVariantReads(variantId, variantContext[0], variantContext[2], variantContext[3], bamFile, True, self.variantContextUnmappedDonor)	# Obtain all donor reads overlapping with the combined variant context and their mates.
								
								# Check whether reads were found in both acceptor and donor. Only then save the results.
								if((len(donorContextReads) > 0) and (len(acceptorContextReads) > 0)):
									variantSampleMap[variantId] = sampleId
									
									self.acceptorContextMap[variantId] = acceptorVariantContext
									self.acceptorContextReadMap[variantId] = acceptorVariantReads
									
									self.donorContextReadMap[variantId] = donorVariantReads
									self.donorContextMap[variantId] = donorVariantContext
									
									self.variantContextMap[variantId] = variantContext
									self.variantContextAcceptorReads[variantId] = acceptorContextReads
									self.variantContextDonorReads[variantId] = donorContextReads
								else:
									self.vaseLogger.debug("No donor and/or acceptor BAM reads found for variant " +str(variantId))
								
							except IOError as ioe:
								self.vaseLogger.warning("Could not obtain BAM reads from " +bamSampleMap[sampleId])
						else:
							self.vaseLogger.debug("VCF variant " +str(variantId)+ " is located in an already existing variant context")
					bamFile.close()
					vcfFile.close()
					
					donorVcfsUsed.append(vcfSampleMap[sampleId])
					donorBamsUsed.append(bamSampleMap[sampleId])
				except IOError as ioe:
					self.vaseLogger.warning("Could not establish data for " +str(sampleId))
			
			
			# Write the variant context data and used donor VCFs/BAMs to output files.
			self.vaseLogger.info("Writing combined variant contexts to " +str(varConOutPath))
			self.writeVariantContexts(self.variantContextMap, self.acceptorContextMap, self.variantContextAcceptorReads, self.donorContextMap, self.variantContextDonorReads, variantSampleMap, varConOutPath)	#Writes the variant contexts to file.
			self.vaseLogger.info("Writing acceptor variant contexts to " +str(varConOutPath)+"/acceptorcontexs.txt")
			self.writeADVariantContext(self.acceptorContextMap, self.acceptorContextReadMap, variantSampleMap, outPath+"/acceptorcontexts.txt")	# Write all acceptor variant contexts to 'acceptorcontexts.txt'
			self.vaseLogger.info("Writing donor variant contexts to " +str(varConOutPath)+"/donorcontexts.txt")
			self.writeADVariantContext(self.donorContextMap, self.donorContextReadMap, variantSampleMap, outPath+"/donorcontexts.txt")	# Write all donor variant contexts to 'donorcontexts.txt'
			
			# Write the read identifiers of reads with unmapped mates to file.
			self.vaseLogger.info("Writing the ids of acceptor reads overlapping with the combined variant context that have unmapped mates to " +str(outPath)+"/varcon_unmapped_acceptor.txt")
			self.writeUnmappedMateReads(self.variantContextUnmappedAcceptor, outPath+"/varcon_unmapped_acceptor.txt")	# All acceptor reads overlapping the combined variant context that have unmapped mates.
			self.vaseLogger.info("Writing the ids of donor reads overlapping with the combined variant context that have unmapped mates to " +str(outPath)+"/varcon_unmapped_donor.txt")
			self.writeUnmappedMateReads(self.variantContextUnmappedDonor, outPath+"/varcon_unmapped_donor.txt")	# All donor reads overlapping the combined variant context that have unmapped mates.
			self.vaseLogger.info("Writing the ids of acceptor reads overlapping with the acceptor variant context to " +str(outPath)+"/acceptor_unmapped.txt")
			self.writeUnmappedMateReads(self.acceptorUnmappedMateMap, outPath+"/acceptor_unmapped.txt")	# Acceptor reads overlapping the variant that have unmapped mates.
			self.vaseLogger.info("Writing the ids of donor reads overlapping with the donor variant context to " +str(outPath)+"/donor_unmapped.txt")
			self.writeUnmappedMateReads(self.donorUnmappedMateMap, outPath+"/donor_unmapped.txt")	# Donor reads overlapping the variant that have unmapped mates.
			
			# Write the new statistics output files.
			self.vaseLogger.info("Writing combined variant context statistics to " +str(outPath)+"/varconstats.txt")
			self.writeVariantContextStats(self.variantContextMap, self.variantContextAcceptorReads, self.variantContextDonorReads, outPath+"/varconstats.txts")	# Writes average and median statistics per variant context
			self.vaseLogger.info("Writing left and right most positions of all acceptor reads and mates per variant context to " +str(outPath)+"/varcon_acceptor_positions.txt")
			self.writeLeftRightPositions(self.variantContextAcceptorReads, outPath+"/varcon_acceptor_positions.txt")	# Write the left (for R1) and right (for R2) most acceptor read positions per variant context to file.
			self.vaseLogger.info("Writing left asnd right most positions of all donor reads and mates per variant context to " +str(outPath)+"/varcon_donor_positions.txt")
			self.writeLeftRightPositions(self.variantContextDonorReads, outPath+"/varcon_donor_positions.txt")	# Write the left (for R1) and right (for R2) most donor read positions per variant context to file.
			
			# Write the used VCF/BAM donor files.
			self.vaseLogger.info("Write the used donor VCF files per sample to " +str(outPath)+"/donorvcfs.txt")
			self.writeUsedDonorFiles(outPath+"/donorvcfs.txt", vcfSampleMap, donorVcfsUsed)
			self.vaseLogger.info("Write the used donor BAM files per sample to " +str(outPath)+"/donorbams.txt")
			self.writeUsedDonorFiles(outPath+"/donorbams.txt", bamSampleMap, donorBamsUsed)
			
			
			# Obtain a list of NIST reads to skip when iterating over the NIST FastQ.
			acceptorReadsToSkip = self.makeTemplateExludeList(self.variantContextAcceptorReads)	# Set up a list of all NIST reads to skip.
			
			# Make the new FastQ files that can be used to run in the NGS_DNA pipeline along real sample data
			self.vaseLogger.info("Start writing the R1 FastQ files")
			self.buildFastQ(fastqFPath, acceptorReadsToSkip, self.donorContextReadMap, 'F', fastqOutPath)	# Build the R1 fastq file.
			self.vaseLogger.info("Wrote all R1 FastQ files")
			
			self.vaseLogger.info("Start writing the R2 FastQ files")
			self.buildFastQ(fastqRPath, acceptorReadsToSkip, self.donorContextReadMap, 'R', fastqOutPath)	# Build the R2 fastq file.
			self.vaseLogger.info("Wrote all R2 FastQ files")
			
			self.vaseLogger.info("Finished building the validation set")
			acceptorBamFile.close()
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not open acceptor BAM file")
			exit()
	
	
	# Returns the BAM reads containing the specific vcf variant as well as their read mate.
	def getVariantReads(self, contextid, vcfVariantChr, varStartPos, varEndPos, bamFile, writeUnm=False, uMateMap=None):
		# Obtain all the variant reads overlapping with the variant and their mate reads.
		variantReads = []
		uMateMap[contextid] = []	#Create an entry for the current contextid
		
		for vread in bamFile.fetch(vcfVariantChr, varStartPos, varEndPos):
			variantReads.append(DonorBamRead(vread.query_name, self.getReadPairNum(vread), vread.reference_name, vread.reference_start, vread.infer_read_length(), vread.get_forward_sequence(), ''.join([chr((x+33)) for x in vread.get_forward_qualities()]), vread.mapping_quality))
			
			# Try to obtain the reads mate as well.
			try:
				vmread = bamFile.mate(vread)
				variantReads.append(DonorBamRead(vmread.query_name, self.getReadPairNum(vmread), vmread.reference_name, vmread.reference_start, vmread.infer_read_length(), vmread.get_forward_sequence(), ''.join([chr((x+33)) for x in vmread.get_forward_qualities()]), vmread.mapping_quality))
			except ValueError as pve:
				self.vaseLogger.debug("Could not find mate for " +vread.query_name+ " ; mate is likely unmapped.")
				if(writeUnm):
					uMateMap[contextid].append(vread.query_name)
		
		# Make sure the list only contains each BAM read once (if a read and mate both overlap with a variant, they have been added twice to the list)
		variantReads = list(set(variantReads))
		variantReads = self.filterVariantReads(variantReads)
		self.vaseLogger.debug("Found a total of " + str(len(variantReads))+ " BAM reads.")
		return variantReads
	
	
	# Filters the donor reads to keep only reads that occur twice.
	def filterVariantReads(self, bamReads):
		filteredList = []
		for bread in bamReads:
			# Add the read to the new list if there are two reads with the same ID (they are a correct pair in that case)
			if(sum(sumread.getBamReadId() == bread.getBamReadId() for sumread in bamReads) == 2):
				filteredList.append(bread)
		return filteredList
	
	
	# Determines the start and stops of the variant context (please see the documentation for more information).
	def determineContext(self, bamVariantReads, contextOrigin):
		# Check whether there are reads to determine the context for.
		if(len(bamVariantReads) > 0):
			# First determine the context start by sorting the reads on leftmost position in ascending order.
			bamVariantReads.sort(key=lambda x:x.getBamReadRefPos(), reverse=False)
			contextChrom = bamVariantReads[0].getBamReadChrom()
			contextStart = bamVariantReads[0].getBamReadRefPos()
			
			# Second determine the context stop by iterating over the reads and calculating the rightmost position of the reads.
			contextStop = 0
			for bvRead in bamVariantReads:
				if(bvRead.getBamReadLength() is not None):
					stopPos = bvRead.getBamReadRefEnd()
					if(stopPos > contextStop):
						contextStop = stopPos
			
			self.vaseLogger.debug("Context is " +str(contextChrom)+ ", " +str(contextStart)+ ", " +str(contextStop))
			return [contextChrom, contextOrigin, contextStart, contextStop]
		return []
	
	
	# Returns whether a certain variant is in an already established variant context.
	def isInContext(self, vcfVarChrom, vcfVarPos):
		for vcfVar, context in self.variantContextMap.items():
			if(vcfVarChrom == context[0]):
				if(vcfVarPos >= context[2] and vcfVarPos <= context[3]):
					return True
		return False
	
	
	# Returns whether a variant has already been used before in the analysis.
	def variantAlreadyProcessed(self, vcfVarId):
		return vcfVarId in self.variantContextMap
	
	
	# Will build the R1/R2 VaSe fastq files.
	def buildFastQ(self, acceptorFqFilePaths, acceptorReadsToSkip, donorContextReadMap, forwardOrReverse, vaseFqOutPath):
		writeDonor = False
		
		# Iterate over the R1/R2 fastq in files to use as templates for the 
		for x in range(0, len(acceptorFqFilePaths)):
			if(x == len(acceptorFqFilePaths)-1):
				writeDonor = True
				self.vaseLogger.debug("Donor reads will be added the current VaSe fastQ out file.")
			
			# Write the new VaSe FastQ file
			vaseFqOutName = self.setFastqOutPath(vaseFqOutPath, forwardOrReverse, x+1)
			self.writeVaSeFastQ(acceptorFqFilePaths[x], vaseFqOutName, acceptorReadsToSkip, donorContextReadMap, forwardOrReverse, writeDonor)
	
	
	# Builds a new FastQ file to be used for validation.
	def writeVaSeFastQ(self, acceptorFastqIn, fastqOutPath, acceptorReadsToSkip, donorBamReadData, fR, writeDonorData=False):
		try:
			fqFile = io.BufferedWriter(gzip.open(fastqOutPath, 'wb'))
			self.vaseLogger.debug("Opened template FastQ: " +acceptorFastqIn)
			
			#Open the template fastq and write filtered data to a new fastq.gz file
			gzFile = io.BufferedReader(gzip.open(acceptorFastqIn, 'rb'))
			for fileLine in gzFile:
				
				# Check if we are located at a read identifier
				if(fileLine.startswith(b"@")):
					if(fileLine.decode("utf-8").strip()[1:] not in acceptorReadsToSkip):
						fqFile.write(fileLine)
						fqFile.write(next(gzFile))
						fqFile.write(next(gzFile))
						fqFile.write(next(gzFile))
			gzFile.close()
			
			# Add the patient BAM reads containing a VCF variant to the new FastQ file.
			if(writeDonorData):
				for vcfvar in donorBamReadData:
					donorBamReads = donorBamReadData[vcfvar]
					donorBamReads.sort(key=lambda x: x.getBamReadId(), reverse=False)
					for bamRead in donorBamReads:
						# Check if the BAM read is R1 or R2.
						if(self.isRequiredRead(bamRead, fR)):
							fqFile.write(bamRead.getAsFastQSeq().encode("utf-8"))
			fqFile.flush()
			fqFile.close()
			
		except IOError as ioe:
			if(ioe.filename==acceptorFastqIn):
				self.vaseLogger.critical("The supplied template FastQ file could not be found.")
			if(ioe.filename==fastqOutPath):
				self.vaseLogger.critical("A FastQ file could not be written to the provided output location.")
			exit()

	
	
	# Checks if a read is read 1 (R1) or read 2 (R2).
	def isRequiredRead(self, bamRead, fR):
		if(fR=="F"):
			if(bamRead.isRead1()):
				return True
		else:
			if(bamRead.isRead2()):
				return True
		return False
	
	
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
	
	
	# Returns the context start and stop for a specified VCF variant.
	def getVariantContext(self, vcfVariant):
		if(vcfVariant in self.variantContextMap):
			return self.variantContextMap[vcfVariant]
		else:
			return None
	
	
	# Adds a VCF variant and its context start and stop if it isn't hasn't been established before.
	def addVariantContext(self, vcfVariant, contextStart, contextStop):
		if(vcfVariant not in self.variantContextMap):
			self.variantContextMap[vcfVariant] = [contextStart, contextStop]
	
	
	#Writes the NIST BAM reads associated with variants to a specified output file.
	def writeAcceptorBamReads(self, nistVariantReadMap, nistBreadOutPath):
		try:
			with open(nistBreadOutPath, 'w') as nistBreadFile:
				nistBreadFile.write("#Variant\tReads\n")
				for variant, acceptorReads in nistVariantReadMap.items():
					nistBreadFile.write(variant + "\t" + " ; ".join([str(x.getBamReadId()) for x in acceptorReads]) + "\n")
				nistBreadFile.close()
		except IOError as ioe:
			self.vaseLogger.critical("Could not write acceptor variant reads to " +ioe.filename)
			exit()
	
	
	# Returns the name for the fastq out file.
	def setFastqOutPath(self, outPath, fR, lNum):
		valName = outPath+ "_" +str(datetime.now().date())+ "_L" +str(lNum)+ "_R2.fastq.gz"
		if(fR=="F"):
			valName = outPath+ "_" +str(datetime.now().date())+ "_L" +str(lNum)+ "_R1.fastq.gz"
		self.vaseLogger.debug("Set FastQ output path to: " +str(valName))
		return valName
	
	
	# Returns an identifier for a VCF variant. If the identifier is '.' then one will be constructed as 'chrom_pos'. 
	def getVcfVariantId(self, vcfVariant):
		return str(vcfVariant.chrom) +"_"+ str(vcfVariant.pos)
	
	
	# Returns the number of occurences of a certain read in the list of BAM reads (should be two ideally)
	def readOccurence(self, readId, readsList):
		return sum(sumread.query_name == readId.query_name for sumread in readsList)
	
	
	# Writes the identifiers of reads that have unmapped mates per sample to a file. Samples are all donors and the ?template?.
	def writeReadsWithUnmappedMates(self, mapOfUnmappedMates, umFileLoc):
		try:
			with open(umFileLoc, 'w') as umFile:
				self.vaseLogger.info("Start writing read ids with unmapped mates to " +umFileLoc)
				umFile.write("#Sample\tRead_IDs\n")
				for sampleId, readIdList in mapOfUnmappedMates.items():
					if(sampleId != "template"):
						umFile.write(sampleId +"\t"+ ";".join(readIdList)+ "\n")
				self.vaseLogger.info("Wrote read ids with unmapped mates to " +umFileLoc)
		except IOError as ioe:
			self.vaseLogger.critical("Could not write ids of reads with unmapped mates to " +umFileLoc)
			exit()
	
	
	# Returns whether a variant is a SNP or indel.
	def determineVariantType(self, vcfVariantRef, vcfVariantAlts):
		maxAltLength = 0
		
		# Determine the maximum length of the alternative allele(s).
		for altAllele in list(vcfVariantAlts):
			if(len(altAllele) > maxAltLength):
				maxAltLength = len(altAllele)
		
		# Check based on the reference and alternative lengths whether the variant is a SNP or indel.
		if(len(vcfVariantRef)==1 and maxAltLength==1):
			return "snp"
		elif(len(vcfVariantRef)>1 or maxAltLength>1):
			return "indel"
		return "?"
	
	
	# Returns the search start and stop to use for searching BAM reads overlapping with the range of the indel
	def determineIndelReadRange(self, variantPos, variantRef, variantAlts):
		searchStart = variantPos
		searchStop = variantPos + len(variantRef)
		
		for varalt in list(variantAlts):
			if((variantPos + len(varalt)) > searchStop):
				searchStop = variantPos + len(varalt)
		return [searchStart, searchStop]
	
	
	# Returns a list of identifiers to skip when making the new validation fastqs
	def makeTemplateExludeList(self, nistReadMap):
		templExcludeList =  []
		for vcVar, templReads in nistReadMap.items():
			for bread in templReads:
				if(bread is not None):
					if(bread.getBamReadId() not in templExcludeList):
						templExcludeList.append(bread.getBamReadId())
		return templExcludeList
	
	
	# Writes the used donor vcf files to a file
	def writeUsedDonorFiles(self, outLocFile, fileSampleMap, listOfUsedDonorFiles):
		try:
			with open(outLocFile, 'w') as outFile:
				outFile.write("#SampleId\tDonorFile")
				for sampleid, sampleFile in fileSampleMap.items():
					if(sampleFile in listOfUsedDonorFiles):
						outFile.write(sampleid+ "\t" +sampleFile+ "\n")
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not write used donor files to " +str(outLocFile))
	
	
	# Writes the read identifiers with unmapped mates to a file.
	def writeUnmappedMateReads(self, unmappedMateMap, outLocFile):
		try:
			with open(outLocFile, 'w') as outFile:
				outFile.write("VariantContext\tReads\n")
				for contextid, ureads in unmappedMateMap.items():
					#unmapreads = self.getVariantContextReadIds(contextid, ureads)
					outFile.write(contextid +"\t"+ ";".join(ureads) +"\n")
		except IOError as ioe:
			self.vaseLogger.warning("Could not write the unmapped mates to file")
	
	
	# Returns whether the read is the first or second read in a pair.
	def getReadPairNum(self, bamRead):
		if(bamRead.is_read1):
			return '1'
		return '2'
	
	
	# Determines the size of the variant context based on both the acceptor and donor reads
	def determineLargestContext(self, contextOrigin, acceptorContext, donorContext):
		largestContext = [acceptorContext[0]]
		largestContext.append(contextOrigin)
		
		# Check the smallest context start
		contextStart = donorContext[2]
		if(acceptorContext[2] < donorContext[2]):
			contextStart = acceptorContext[2]
		largestContext.append(contextStart)
		
		# Check the largest context end
		contextEnd = donorContext[3]
		if(acceptorContext[3] > donorContext[3]):
			contextEnd = acceptorContext[3]
		largestContext.append(contextEnd)
		return largestContext
	
	
	
	# ====================NEWLY ADDED METHODS====================
	
	# Determines the start and end positions to use for searching reads overlapping with the variant
	def determineReadSearchWindow(self, variantType, vcfVariant):
		if(variantType == 'snp'):
			return [vcfVariant.pos-1, vcfVariant.pos+1]
		elif(variantType == 'indel'):
			return self.determineIndelReadRange(vcfVariant.ref, vcfVariant.alts)
		return [-1, -1]
	
	
	# Writes the new variant context output file format that combines the old varcon with the acceptorbread and donorbread files into one
	def writeVariantContexts(self, variantContextMap, acceptorContextMap, acceptorContextReads, donorContextMap, donorContextReads, variantSampleMap, varConOutPath):
		try:
			with open(varConOutPath, 'w') as varconFile:
				varconFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tAcceptorContext\tDonorContext\tAcceptorReads\tDonorReads\tADratio\tAcceptorReadsIds\tDonorReadIds\n")
				
				# Iterate over all the variants and their contexts.
				for variant, varContext in variantContextMap.items():
					varconAreads = self.getVariantContextReadIds(variant, acceptorContextReads)
					varconDreads = self.getVariantContextReadIds(variant, donorContextReads)
					adRatio = len(varconAreads) / len(varconDreads)
					varconFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+  varContext[0] +"\t"+ str(varContext[1]) +"\t"+ str(varContext[2]) +"\t"+ str(varContext[3]) +"\t"+ str(len(acceptorContextMap[variant])) +"\t"+ str(len(donorContextMap[variant])) +"\t"+ str(len(varconAreads)) +"\t"+ str(len(varconDreads)) +"\t"+ str(adRatio) +"\t"+ ';'.join(varconAreads) +"\t"+ ';'.join(varconDreads) +"\n")
			self.vaseLogger.info("Finished writing variants and their contexts to " +varConOutPath)
		except IOError as ioe:
			self.vaseLogger.critical("Could not write variant contexts to " +str(varConOutPath))
			exit()
	
	
	# Writes the variant contexts for acceptor and donor variant contexts.
	def writeADVariantContext(self, adVarconMap, adVarconReads, variantSampleMap, adVarconOutPath):
		try:
			with open(adVarconOutPath, 'w') as advcFile:
				advcFile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tNumOfReads\tReadIds\n")
				
				# Iterate over the acceptor/donor variant contexts.
				for contextId, variantcontext in adVarconMap.items():
					adreads = self.getVariantContextReadIds(contextId, adVarconReads)
					advcFile.write(str(contextId) +"\t"+ str(variantSampleMap[contextId]) +"\t"+ str(variantcontext[0]) +"\t"+ str(variantcontext[1]) +"\t"+ str(variantcontext[2]) +"\t"+ str(variantcontext[3]) +"\t"+ str(len(adreads)) +"\t"+ ';'.join(adreads))
			self.vaseLogger.info("Finished writing acceptor/donor variant contexts to " +str(adVarconOutPath))
		except IOError as ioe:
			self.vaseLogger.critical("Could not write acceptor/donor variant contexts to " +str(adVarconOutPath))
	
	
	# Returns the saved read objects for a specific variant context
	def getVariantContextReads(self, contextId, contextReadMap):
		if(contextId in contextReadMap):
			return contextReadMap[contextId]
		return []
	
	
	# Returns the list of read identifiers for a specified variant context and read map.
	def getVariantContextReadIds(self, contextId, contextReadMap):
		readIdList = []
		if(contextId in contextReadMap):
			for readObj in contextReadMap[contextId]:
				readIdList.append(readObj.getBamReadId())
		return readIdList
	
	
	# Writes the left most positions (to be able to construct a histogram)
	def writeLeftMostPositionPerVariantContext(self, varconReadData, leftPosOutFileLoc):
		try:
			with open(leftPosOutFileLoc, 'w') as lpof:
				lpof.write("#ContextId\tLeftPos\n")
				for contextId, varconReads in varconReadData.items():
					leftPositions = []
					for vcRead in varconReads:
						if(vcRead.getBamReadRefPos() is None):
							print(str(vcRead.getBamReadId) +" ; "+ str(vcRead.getBamReadPairNumber()))
						leftPositions.append(str(vcRead.getBamReadRefPos()))
					lpof.write(str(contextId)+ "\t" +','.join(leftPositions)+ "\n")
		except IOError as ioe:
			self.vaseLogger.warning("Could not write read left positions to output file " +str(leftPosOutFileLoc))
	
	
	# Writes the left and right positions to the output file. Left pos for R1 and right pos for R2
	def writeLeftRightPositions(self, varconReadData, outFileLoc):
		try:
			with open(outFileLoc, 'r') as lrpof:
				lrpof.write("#ContextId\tLeftPos\tRightPos\n")
				for contextid, varconReads in varconReadData.items():
					leftPositions, rightPositions = [], []
					
					# Write the left and right positions for each variant context id
					for vcRead in varconReads:
						if(vcRead.isRead1()):
							leftPositions.append(str(vcRead.getBamReadRefPos()))
						else:
							rightPositions.append(str(vcRead.getBamReadRefEnd()))
					lrpof.write(str(contextId) +"\t"+ ','.join(leftPositions) +"\t"+ ','.join(rightPositions))
		except IOError as ioe:
			self.vaseLogger.warning("Could not write read left positions to output file " +str(outFileLoc))
	
	
	# Writes some statistics about the acceptor and donor reads identified for each variant context.
	def writeVariantContextStats(self, variantcontexts, varconacceptorreads, varcondonorreads, outFileLoc):
		try:
			with open(outFileLoc, 'w') as statsOutFile:
				statsOutFile.write("#ContextId\tAVG_ALength\tAVG_DLength\tMEDIAN_ALength\tMEDIAN_DLength\tAVG_AQual\tAVG_DQual\tMEDIAN_AQual\tMEDIAN_DQual\tAVG_AMapQ\tAVG_DMapQ\tMEDIAN_MapQ\tMEDIAN_DMapQ\n")
				for contextid in variantcontexts:
					acceptorLengthData = self.getAverageAndMedianReadLength(contextid, varconacceptorreads)
					acceptorQualData = self.getAverageAndMedianReadQual(contextid, varconacceptorreads)
					acceptorMapQData = self.getAverageAndMedianReadMapQ(contextid, varconacceptorreads)
					donorLengthData = self.getAverageAndMedianReadLength(contextid, varcondonorreads)
					donorQualData = self.getAverageAndMedianReadQual(contextid, varcondonorreads)
					donorMapQData = self.getAverageAndMedianReadMapQ(contextid, varcondonorreads)
					
					# Write the data for the variant context to file.
					statsOutFile.write(str(contextid)+ "\t" +str(acceptorLengthData[0])+ "\t" +str(donorLengthData[0])+ "\t" +str(acceptorLengthData[1])+ "\t" +str(donorLengthData[1])+ "\t" +str(acceptorQualData[0])+ "\t" +str(donorQualData[0])+ "\t" +str(acceptorQualData[1])+ "\t" +str(donorQualData[1])+ "\t" +str(acceptorMapQData[0])+ "\t" +str(donorMapQData[0])+ "\t" +str(acceptorMapQData[1])+ "\t" +str(donorMapQData[1])+ "\n")
		except IOError as ioe:
			self.vaseLogger.warning("Could not write variant context statistics to " +str(outFileLoc))
	
	
	# Returns the average and median read quality for a variant context
	def getAverageAndMedianReadQual(self, contextid, variantcontextreads):
		if(contextid in variantcontextreads):
			avgMedQual = []
			for contextread in variantcontextreads[contextid]:
				avgMedQual.append(contextread.getAverageQscore())
			return ([statistics.mean(avgMedQual), statistics.median(avgMedQual)])
		return None
	
	
	# Returns the average and median read mapq value for a variant context
	def getAverageAndMedianReadMapQ(self, contextid, variantcontextreads):
		if(contextid in variantcontextreads):
			avgMedMapQ = []
			for contextread in variantcontextreads[contextid]:
				avgMedMapQ.append(contextread.getMappingQual())
			return ([statistics.mean(avgMedMapQ), statistics.median(avgMedMapQ)])
		return None
	
	
	# Returns the average and median read length for a variant context
	def getAverageAndMedianReadLength(self, contextid, variantcontextreads):
		if(contextid in variantcontextreads):
			avgMedLen = []
			for contextread in variantcontextreads[contextid]:
				avgMedLen.append(contextread.getBamReadLength())
			return ([statistics.mean(avgMedLen), statistics.median(avgMedLen)])
		return None
