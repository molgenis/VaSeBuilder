#!/usr/bin/env python
import io
import logging
from datetime import datetime
import gzip
import pysam
from DonorBamRead import DonorBamRead

class VaSeBuilder:
	# Constructor that saves the identifier, date and time of the current 
	def __init__(self, vaseId):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.creationId = str(vaseId)
		self.creationDate = datetime.now().date()
		self.creationTime = datetime.now().time()
		self.variantContextMap = {}	# Saves the variant contexts as {variantId => [chrom, start, stop]}
		self.variantBamReadMap = {}	# Saves the donor BAM reads per variant context as {variantId => [dbamReadId1, dbamReadId2, ..., dbamReadIdn]}
		self.nistVariantReadMap = {}	# Saves the acceptor BAM reads per variant context as {variantId => [abamReadId1, abamReadId2, ..., abamReadIdn]}
		self.variantBamFileMap = {}	# Saves which BAM file to use per variant (not needed anymore! but too lazy to remove)
		self.unmappedMateMap = {}	# Saves the BAM read identifiers that have an unmapped mate per sample as {sampleId => [ubamRead1, ubamRead2, ..., ubamReadn]}
		self.vaseLogger.info("VaSeBuilder: " +self.creationId+ " ; " +str(self.creationDate)+ " ; " +str(self.creationTime))
	
	
	# Creates the new FastQ validation dataset by replacing NIST reads containing a VCF variant with BAM reads from patients. Returns true at the end to indicate the process is done.
	def buildValidationSet(self, vcfBamLinkMap, vcfSampleMap, bamSampleMap, nistBamLoc, fastqFPath, fastqRPath, outPath, fastqOutPath, varConOutPath, varBreadOutPath, nistBreadOutPath):
		self.vaseLogger.info("Start building the validation set")
		variantSampleMap = {}	# Will link the variants to their samples. (Is used when writing the variant data to files)
		donorVcfsUsed, donorBamsUsed = [], []
		self.unmappedMateMap['template'] = []
		
		try:
			acceptorBamFile = pysam.AlignmentFile(nistBamLoc, 'rb')
			
			# Iterate over the samples to use for building the validation set.
			for sampleId in vcfSampleMap:
				self.vaseLogger.debug("Processing data for sample " +sampleId)
				self.unmappedMateMap[sampleId] = []
				
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
						variantSampleMap[variantId] = sampleId
						
						# Get the BAM reads fr the variant and determine the variant context.
						if(not self.isInContext(vcfVar.chrom, vcfVar.pos)):
							try:
								variantReads, variantContext, acceptorReads = [], [], []
								searchStart, searchStop = 0, 0
								
								# Determine the variant type before gathering reads overlapping with the VCF variant
								variantType = self.determineVariantType(vcfVar.ref, vcfVar.alts)
								if(variantType == 'snp'):
									searchStart, searchStop = vcfVar.pos-1, vcfVar.pos+1
								elif(variantType == 'indel'):
									indelPositions = self.determineIndelReadRange(vcfVar.ref, vcfVar.alts)
									searchStart, searchStop = indelPositions[0], indelPositions[1]
								
								self.vaseLogger.debug("Search donor BAM reads for variant " +variantId)
								variantReads = self.getVariantReads(sampleId, vcfVar.chrom, searchStart, searchStop, bamFile)	# Obtain all patient BAM reads containing the VCF variant and their read mate.
								
								self.vaseLogger.debug("Determine context for variant " +variantId)
								variantContext = self.determineContext(variantReads)	# Save the context start and stop for the variant in the VCF variant context map.
								variantContext.insert(1, vcfVar.pos)	# Add the position of the variant the variant context is based on
								
								self.vaseLogger.debug("Determine template BAM reads for variant " +variantId)
								acceptorReads = self.getVariantReads("template", vcfVar.chrom, searchStart, searchStop, acceptorBamFile)	# Obtain all NIST BAM reads containing the VCF variant. and their read mate.
								
								
								# Check whether reads were found in both patient and NIST. Only then save the results.
								if((len(variantReads) > 0) and (len(acceptorReads) > 0)):
									self.variantBamReadMap[variantId] = variantReads
									self.variantContextMap[variantId] = variantContext
									self.nistVariantReadMap[variantId] = acceptorReads
									self.variantBamFileMap[variantId] = bamSampleMap[sampleId]
								else:
									self.vaseLogger.debug("No donor or template BAM reads found for variant " +variantId)
								
							except IOError as ioe:
								self.vaseLogger.warning("Could not obtain BAM reads from " +bamSampleMap[sampleId])
						else:
							self.vaseLogger.debug("VCF variant " +variantId+ " is located in an already existing variant context")
					bamFile.close()
					vcfFile.close()
					
					donorVcfsUsed.append(vcfSampleMap[sampleId])
					donorBamsUsed.append(bamSampleMap[sampleId])
				except IOError as ioe:
					self.vaseLogger.warning("Could not establish data for " +str(sampleId))
			
			# Write data used to build the new FastQ to output files.
			self.writeVariantsContexts(self.variantContextMap, variantSampleMap, varConOutPath)	# Write the context start and stop for each used variant to a separate file.
			self.writeVariantBamReads(self.variantBamFileMap, self.variantBamReadMap, variantSampleMap, varBreadOutPath)	# Write the associated BAM reads for each used variant to a seperate file.
			#self.writeNistVariantBamReads(self.nistVariantReadMap, nistBreadOutPath)	# Write the associated NIST BAM reads for each used variant to a separate file.
			self.writeAcceptorBamReads(self.nistVariantReadMap, nistBreadOutPath)	# Write the associated NIST BAM reads for each used variant to a separate file.
			
			# Obtain a list of NIST reads to skip when iterating over the NIST FastQ.
			acceptorReadsToSkip = self.makeTemplateExludeList(self.nistVariantReadMap)	# Set up a list of all NIST reads to skip.
			
			# Write all used donor VCF and BAM files to output files.
			self.writeUsedDonorFiles(outPath+"/donorvcfs.txt", vcfSampleMap, donorVcfsUsed)
			self.writeUsedDonorFiles(outPath+"/donorbams.txt", bamSampleMap, donorBamsUsed)
			self.writeUnmappedMateReads(self.unmappedMateMap, outPath+"/unmappedmatereads.txt")
			
			
			# Make the new FastQ files that can be used to run in the NGS_DNA pipeline along real sample data
			self.vaseLogger.info("Start writing the first (_R1) FastQ file")
			self.buildFastQ(fastqFPath, self.setFastqOutPath(fastqOutPath, 'F'), acceptorReadsToSkip, self.variantBamReadMap, self.variantBamFileMap, 'F')	# Build the R1 fastq file.
			self.vaseLogger.info("Wrote the first (_R1) FastQ file")
			
			self.vaseLogger.info("Start writing the second (_R2) FastQ file")
			self.buildFastQ(fastqRPath, self.setFastqOutPath(fastqOutPath, 'R'), acceptorReadsToSkip, self.variantBamReadMap, self.variantBamFileMap, 'R')	# Build the R2 fastq file.
			self.vaseLogger.info("Wrote the second (_R2) FastQ file")
			
			self.vaseLogger.info("Finished building the validation set")
			acceptorBamFile.close()
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not open acceptor BAM file")
			exit()
	
	
	# Returns the BAM reads containing the specific vcf variant as well as their read mate.
	def getVariantReads(self, sampleid, vcfVariantChr, varStartPos, varEndPos, bamFile):
		# Obtain all the variant reads overlapping with the variant and their mate reads.
		variantReads = []
		for vread in bamFile.fetch(vcfVariantChr, varStartPos, varEndPos):
			variantReads.append(DonorBamRead(vread.query_name, self.getReadPairNum(vread), vread.reference_name, vread.reference_start, vread.infer_read_length(), vread.get_forward_sequence(), ''.join([chr((x+33)) for x in vread.get_forward_qualities()])))
			
			# Try to obtain the reads mate as well.
			try:
				vmread = bamFile.mate(vread)
				variantReads.append(DonorBamRead(vmread.query_name, self.getReadPairNum(vmread), vmread.reference_name, vmread.reference_start, vmread.infer_read_length(), vmread.get_forward_sequence(), ''.join([chr((x+33)) for x in vmread.get_forward_qualities()])))
			except ValueError as pve:
				self.vaseLogger.debug("Could not find mate for " +vread.query_name+ " ; mate is likely unmapped.")
				self.unmappedMateMap[sampleid].append(vread.query_name)
		
		# Make sure the list only contains each BAM read once (if a read and mate both overlap with a variant, they have been added twice to the list)
		variantReads = list(set(variantReads))
		variantReads = self.filterDonorReads(variantReads)
		self.vaseLogger.debug("Found a total of " + str(len(variantReads))+ " BAM reads.")
		return variantReads
	
	
	# Filters the donor reads to keep only reads that occur twice.
	def filterDonorReads(self, bamReads):
		filteredList = []
		for bread in bamReads:
			# Add the read to the new list if there are two reads with the same ID (they are a correct pair in that case)
			if(sum(sumread.getBamReadId() == bread.getBamReadId() for sumread in bamReads) == 2):
				filteredList.append(bread)
		return filteredList
	
	
	# Determines the start and stops of the variant context (please see the documentation for more information).
	def determineContext(self, bamVariantReads):
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
			return [contextChrom, contextStart, contextStop]
		return []
	
	
	# Returns whether a certain variant is in an already established variant context.
	def isInContext(self, vcfVarChrom, vcfVarPos):
		for vcfVar, context in self.variantContextMap.items():
			if(vcfVarChrom == context[0]):
				if(vcfVarPos >= context[1] and vcfVarPos <= context[2]):
					return True
		return False
	
	
	# Returns whether a variant has already been used before in the analysis.
	def variantAlreadyProcessed(self, vcfVarId):
		return vcfVarId in self.variantContextMap
	
	
	# Builds the new FastQ file to be used.
	def buildFastQ(self, nistFastqIn, fastqOutPath, acceptorReadsToSkip, patientBamReads, variantBamFileMap, fR):
		try:
			fqFile = io.BufferedWriter(gzip.open(fastqOutPath, 'wb'))
			self.vaseLogger.debug("Opened template FastQ: " +nistFastqIn)
			
			#Open the template fastq and write filtered data to a new fastq.gz file
			gzFile = io.BufferedReader(gzip.open(nistFastqIn, 'rb'))
			for fileLine in gzFile:
				
				# Check if we are located at a read identifier
				if(fileLine.startswith(b"@")):
					if(fileLine.decode("utf-8").strip() not in acceptorReadsToSkip):
						fqFile.write(fileLine)
						fqFile.write(next(gzFile))
						fqFile.write(next(gzFile))
						fqFile.write(next(gzFile))
			gzFile.close()
			
			# Add the patient BAM reads containing a VCF variant to the new FastQ file.
			for vcfvar in patientBamReads:
				donorBamReads = patientBamReads[vcfvar]
				donorBamReads.sort(key=lambda x: x.getBamReadId(), reverse=False)
				for bamRead in donorBamReads:
					# Check if the BAM read is R1 or R2.
					if(self.isRequiredRead(bamRead, fR)):
						fqFile.write(bamRead.getAsFastQSeq().encode("utf-8"))
			fqFile.flush()
			fqFile.close()
			
		except IOError as ioe:
			if(ioe.filename==nistFastqIn):
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
	
	
	# Obtains all required info from the 
	def getBamReadAsFastQ(self, bamRead):
		readName = self.getBamReadFastQName(bamRead)	# Check whether to add /1 or /2 to the read identifier.
		symbolQualities = ''.join([chr((x+33)) for x in bamRead.get_forward_qualities()])	# Convert the Q-Score qualities to ASCII symbols as in the original FastQ file.
		fqEntry = readName + "\n" + bamRead.get_forward_sequence() + "\n+\n" + symbolQualities
		return fqEntry
	
	
	# Returns a proper fastq identifier for a read.
	def getBamReadFastQName(self, bamRead):
		if(not (bamRead.query_name.endswith("/1") or bamRead.query_name.endswith("/2"))):
			if(bamRead.is_read1):
				return "@"+bamRead.query_name+"/1"
			return "@"+bamRead.query_name+"/2"
		return bamRead.query_name
	
	
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
	
	
	# Writes the VCF variants and their contexts to a separate file.
	def writeVariantsContexts(self, variantContextMap, variantSampleMap, varConOutPath):
		try:
			# Write the variants with their contexts to a specified output file.
			self.vaseLogger.info("Start writing variants and their contexts to " +varConOutPath)
			with open(varConOutPath, 'w') as varcoFile:
				varcoFile.write("Variant\tSample\tChrom\tOrigin\tStart\tStop\n")
				
				# Iterate over all the variants and their contexts.
				for variant, varContext in variantContextMap.items():
					varcoFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+  varContext[0] +"\t"+ str(varContext[1]) +"\t"+ str(varContext[2]) +"\n")
			self.vaseLogger.info("Finished writing variants and their contexts to " +varConOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write to " +varConOutPath)
			exit()
	
	
	# Writes the VCF variants and their associated BAM reads to a separate file.
	def writeVariantBamReads(self, variantBamFileMap, variantBamReadMap, variantSampleMap, varBreadOutPath):
		try:
			# Write the variants with their contexts to a specified output file.
			self.vaseLogger.info("Start writing variants and their associated BAM reads to " +varBreadOutPath)
			with open(varBreadOutPath, 'w') as varBreadFile:
				varBreadFile.write("Variant\tSample\tReads\n")
				
				# Iterate over all the variants and their contexts.
				for variant, bamReads in variantBamReadMap.items():
					bamReadIdList = list(set([bread.getBamReadId() for bread in bamReads]))
					varBreadFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+ " ; ".join(bamReadIdList) +"\n")
			self.vaseLogger.info("Finished writing variants and their associated BAM reads to " +varBreadOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write BAM reads to " +varBreadOutPath)
			exit()
	
	
	#Writes the NIST BAM reads associated with variants to a specified output file.
	def writeAcceptorBamReads(self, nistVariantReadMap, nistBreadOutPath):
		try:
			with open(nistBreadOutPath, 'w') as nistBreadFile:
				nistBreadFile.write("Variant\tReads\n")
				for variant, acceptorReads in nistVariantReadMap.items():
					nistBreadFile.write(variant + "\t" + " ; ".join([str(x.getBamReadId()) for x in acceptorReads]) + "\n")
				nistBreadFile.close()
		except IOError as ioe:
			self.vaseLogger.critical("Could not write acceptor variant reads to " +ioe.filename)
			exit()
	
	
	# Returns the name for the fastq out file.
	def setFastqOutPath(self, outPath, fR):
		valName = outPath+ "_" +str(datetime.now().date())+ "_R2.fastq.gz"
		if(fR=="F"):
			valName = outPath+ "_" +str(datetime.now().date())+ "_R1.fastq.gz"
		self.vaseLogger.debug("Set FastQ output path to: " +valName)
		return valName
	
	
	# Returns an identifier for a VCF variant. If the identifier is '.' then one will be constructed as 'chrom_pos'. 
	def getVcfVariantId(self, vcfVariant):
		variantId = vcfVariant.id
		if(variantId=='.' or variantId==None):
			variantId = "SNP"+ str(vcfVariant.chrom) +"_"+ str(vcfVariant.pos)
		return variantId
	
	
	# Adjusts the read info so the read identifier will have either /1 or /2 (SeqIO does not do this even if the read file does have /1 or /2)
	def getSeqIoFastqRead(self, seqioRead, fR):
		seqioReadId = "@"+seqioRead.id
		if(not (seqioRead.id.endswith("/1") or seqioRead.id.endswith("/2"))):
			seqioReadId = "@" +seqioRead.id+ "/2"
			if(fR=="F"):
				seqioReadId = "@" +seqioRead.id+ "/1"
		fqRead = seqioReadId+ "\n" + "\n".join(seqioRead.format("fastq").split("\n")[1:])	# Get the sequence, '+' and quality scores but not the read identifier (this we add ourselves)
		return fqRead
	
	
	# Returns the number of occurences of a certain read in the list of BAM reads (should be two ideally)
	def readOccurence(self, readId, readsList):
		return sum(sumread.query_name == readId.query_name for sumread in readsList)
	
	
	# Writes the identifiers of reads that have unmapped mates per sample to a file. Samples are all donors and the ?template?.
	def writeReadsWithUnmappedMates(self, mapOfUnmappedMates, umFileLoc):
		try:
			with open(umFileLoc, 'w') as umFile:
				self.vaseLogger.info("Start writing read ids with unmapped mates to " +umFileLoc)
				umFile.write("Sample\tRead_IDs\n")
				for sampleId, readIdList in mapOfUnmappedMates.items():
					if(sampleId != "template"):
						umFile.write(sampleId +"\t"+ ";".join(readIdList)+ "\n")
				self.vaseLogger.info("Wrote read ids with unmapped mates to " +umFileLoc)
		except IOError as ioe:
			self.vaseLogger.critical("Could not write ids of reads with unmapped mates to " +umFileLoc)
			exit()
	
	
	# Returns the map containing the identifiers of reads with unmapped mates per sample
	def getUnmappedMateMap(self):
		return self.unmappedMateMap
	
	
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
				for sampleid, sampleFile in fileSampleMap.items():
					if(sampleFile in listOfUsedDonorFiles):
						outFile.write(sampleid+ "\t" +sampleFile+ "\n")
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not write used donor files to " +str(outLocFile))
	
	
	# Writes the read identifiers with unmapped mates to a file.
	def writeUnmappedMateReads(self, unmappedMateMap, outLocFile):
		try:
			with open(outLocFile, 'w') as outFile:
				outFile.write("Sample\tReads\n")
				for sampleid, ureads in unmappedMateMap.items():
					outFile.write(sampleid +"\t"+ " ; ".join(ureads) +"\n")
		except IOError as ioe:
			self.vaseLogger.warning("Could not write the unmapped mates to file")
	
	
	# Returns whether the read is the first or second read in a pair.
	def getReadPairNum(self, bamRead):
		if(bamRead.is_read1):
			return '1'
		return '2'
