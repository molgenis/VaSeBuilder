#!/usr/bin/env python
import logging
from datetime import datetime
from Bio import SeqIO
import gzip
import pysam
import numpy


class VaSeBuilder:
	# Constructor that saves the identifier, date and time of the current 
	def __init__(self, vaseId):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.creationId = vaseId
		self.creationDate = datetime.now().date()
		self.creationTime = datetime.now().time()
		self.variantContextMap = {}
		self.variantBamReadMap = {}
		self.nistVariantReadMap = {}
		self.variantBamFileMap = {}
		self.unmappedMateMap = {}
		self.vaseLogger.info("VaSeBuilder: " +self.creationId+ " ; " +str(self.creationDate)+ " ; " +str(self.creationTime))
	
	
	
	# Creates the new FastQ validation dataset by replacing NIST reads containing a VCF variant with BAM reads from patients. Returns true at the end to indicate the process is done.
	def buildValidationSet(self, vcfSampleMap, bamSampleMap, nistBamLoc, fastqFPath, fastqRPath, fastqOutPath, varConOutPath, varBreadOutPath, nistBreadOutPath):
		self.vaseLogger.info("Start building the valdation set")
		variantSampleMap = {}	# Will link the variants to their samples. (Is used when writing the variant data to files)
		
		try:
			nistBamFile = pysam.AlignmentFile(nistBamLoc, 'rb')
			
			# Iterate over the samples to use for building the validation set.
			for sampleId in vcfSampleMap:
				self.vaseLogger.debug("Processing data for sample " +sampleId)
				
				# Check in the VCF BAM link map that there is a BAM file for the sample as well.
				try:
					vcfFile = pysam.VariantFile(vcfSampleMap[sampleId], 'r')
					self.vaseLogger.debug("Opened VCF file " +vcfSampleMap[sampleId]+ ".")
					
					# Loop over the variants in the VCF file. Prior to identifying the BAM reads, it is first checked whether the variant is in a previously established 
					for vcfVar in vcfFile.fetch():
						variantId = self.getVcfVariantId(vcfVar)
						self.vaseLogger.debug("Searching BAM reads for variant " + variantId)
						variantSampleMap[sampleId] = variantId
						
						# Get the BAM reads fr the variant and determine the variant context.
						if(not self.isInContext(vcfVar.chrom, vcfVar.pos)):
							try:
								bamFile = pysam.AlignmentFile(bamSampleMap[sampleId], 'rb')
								self.vaseLogger.debug("Opened BAM file " +bamSampleMap[sampleId])
								
								self.vaseLogger.debug("Search donor BAM reads for variant " +variantId)
								variantReads = self.getVariantReads(sampleId, vcfVar.chrom, vcfVar.pos, bamFile)	# Obtain all patient BAM reads containing the VCF variant and their read mate.
								
								self.vaseLogger.debug("Determine context for variant " +variantId)
								variantContext = self.determineContext(variantReads)	# Save the context start and stop for the variant in the VCF variant context map.
								
								self.vaseLogger.debug("Determine template BAM reads for variant " +variantId)
								nistReads = self.getVariantReads("template", vcfVar.chrom, vcfVar.pos, nistBamFile)	# Obtain all NIST BAM reads containing the VCF variant. and their read mate.
								
								
								# Check whether reads were found in both patient and NIST. Only then save the results.
								if((len(variantReads) > 0) and (len(nistReads) > 0)):
									self.variantBamReadMap[variantId] = variantReads
									self.variantContextMap[variantId] = variantContext
									self.nistVariantReadMap[variantId] = nistReads
									self.variantBamFileMap[variantId] = bamSampleMap[sampleId]
								else:
									self.vaseLogger.debug("No donor or template BAM reads found for variant " +variantId)
								bamFile.close()
							except IOError as ioe:
								self.vaseLogger.warning("Could not obtain BAM reads from " +bamSampleMap[sampleId])
						else:
							self.vaseLogger.debug("VCF variant " +variantId+ " is located in an already existing variant context")
					vcfFile.close()
				
				except IOError as ioe:
					self.vaseLogger.warning("Could not establish data for " +sampleId)
			
			# Write data used to build the new FastQ to output files.
			self.writeVariantsContexts(self.variantContextMap, variantSampleMap, varConOutPath)	# Write the context start and stop for each used variant to a separate file.
			self.writeVariantBamReads(self.variantBamFileMap, self.variantBamReadMap, variantSampleMap, varBreadOutPath)	# Write the associated BAM reads for each used variant to a seperate file.
			self.writeNistVariantBamReads(self.nistVariantReadMap, nistBreadOutPath)	# Write the associated NIST BAM reads for each used variant to a separate file.
			
			# Obtain a list of NIST reads to skip when iterating over the NIST FastQ.
			nistList = numpy.array(list(self.nistVariantReadMap.values())).ravel()	# Set up a list of all NIST reads to skip.
			nistReadsToSkip = [x.query_name for x in nistList]
			
			
			# Make the new FastQ files that can be used to run in the NGS_DNA pipeline along real sample data
			self.vaseLogger.debug("Start writing the first (_R1) FastQ file")
			self.buildFastQ(fastqFPath, self.setFastqOutPath(fastqOutPath, 'F'), nistReadsToSkip, self.variantBamReadMap, self.variantBamFileMap, 'F')	# Build the R1 fastq file.
			self.vaseLogger.debug("Wrote the first (_R1) FastQ file")
			
			self.vaseLogger.debug("Start writing the second (_R2) FastQ file")
			self.buildFastQ(fastqRPath, self.setFastqOutPath(fastqOutPath, 'R'), nistReadsToSkip, self.variantBamReadMap, self.variantBamFileMap, 'R')	# Build the R2 fastq file.
			self.vaseLogger.debug("Wrote the second (_R2) FastQ file")
			
			self.vaseLogger.info("Finished building the validation set")
			nistBamFile.close()
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not open NIST BAM file")
			exit()
	
	
	
	# Returns the BAM reads containing the specific vcf variant as well as their read mate.
	def getVariantReads(self, sampleid, vcfVariantChr, vcfVariantPos, bamFile):
		# Obtain all the variant reads overlapping with the variant and their mate reads.
		variantReads = []
		for vread in bamFile.fetch(vcfVariantChr, vcfVariantPos-1, vcfVariantPos+1):
			variantReads.append(vread)	# Add the BAM read to the list of reads.
			try:
				variantReads.append(bamFile.mate(vread))	# Add the mate of the current BAM read to the list as well
			except ValueError as pve:
				self.vaseLogger.debug("Could not find mate for " +vread.query_name+ " ; mate is likely unmapped.")
				self.unmappedMateMap[sampleid].append(vread_query_name)
		variantReads = list(set(variantReads))	# Make sure the list only contains each BAM read once (if a read and mate both overlap with a variant, they have been added twice to the list)
		self.vaseLogger.debug("Found a total of " + str(len(variantReads))+ " BAM reads.")
		
		return variantReads
	
	
	
	# Determines the start and stops of the variant context (please see the documentation for more information).
	def determineContext(self, bamVariantReads):
		# Check whether there are reads to determine the context for.
		if(len(bamVariantReads) > 0):
			# First determine the context start by sorting the reads on leftmost position in ascending order.
			bamVariantReads.sort(key=lambda x:x.reference_start, reverse=False)
			contextChrom = bamVariantReads[0].reference_name
			contextStart = bamVariantReads[0].reference_start
			
			# Second determine the context stop by iterating over the reads and calculating the rightmost position of the reads.
			contextStop = 0
			for bvRead in bamVariantReads:
				if(bvRead.infer_read_length() is not None):
					stopPos = bvRead.reference_start + bvRead.infer_read_length()
					if(stopPos > contextStop):
						contextStop = stopPos
			
			self.vaseLogger.debug("Context is " +str(contextChrom)+ ", " +str(contextStart)+ ", " +str(contextStop))
			return [contextChrom, contextStart, contextStop]
		else:
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
	def buildFastQ(self, nistFastqIn, fastqOutPath, nistReadsToSkip, patientBamReads, variantBamFileMap, fR):
		# Build the new NIST FastQ file replacing some of its own reads with patient reads containing VCF variants.
		try:
			# Write NIST fastq entries to the new FastQ file but exclude NIST reads overlapping with the position of a VCF variant.
			with gzip.open(nistFastqIn, "rt") as nistFq, gzip.open(fastqOutPath, 'wt') as newNistFq:
				self.vaseLogger.debug("Opened template FastQ: " +nistFastqIn)
				for fqRead in SeqIO.parse(nistFq, "fastq"):
					if(fqRead.id not in nistReadsToSkip):
						newNistFq.write(self.getSeqIoFastqRead(fqRead, fR))
						#newNistFq.write(fqRead.format("fastq"))
				
				# Add the patient BAM reads containing a VCF variant to the new FastQ file.
				for vcfvar in patientBamReads:
					try:
						bamFile = pysam.AlignmentFile(variantBamFileMap[vcfvar], 'rb')
						donorBamReads = patientBamReads[vcfvar]
						donorBamReads.sort(key=lambda x: x.query_name, reverse=False)
						for bamRead in donorBamReads:
							# Check if the BAM read is R1 or R2.
							if(self.isRequiredRead(bamRead, fR) and self.readOccurence(bamRead, donorBamReads)==2):
								newNistFq.write(self.getBamReadAsFastQ(bamRead)+"\n")
						bamFile.close()
					except IOError as ioe:
						self.vaseLogger.warning("Could not open BAM file to obtain bam read data.")
		
		# Throw an exception if one of the files does not exist.
		except IOError as ioe:
			if(ioe.filename==nistFastqIn):
				self.vaseLogger.critical("The supplied NIST FastQ file could not be found.")
			if(ioe.filename==fastqOutPath):
				self.vaseLogger.critical("A FastQ file could not be written to the provided output location.")
			exit()
	
	
	
	# Checks if a read is read 1 (R1) or read 2 (R2).
	def isRequiredRead(self, bamRead, fR):
		if(fR=="F"):
			if(bamRead.is_read1):
				return True
		else:
			if(bamRead.is_read2):
				return True
		return False
	
	
	
	# Obtains all required info from the 
	def getBamReadAsFastQ(self, bamRead):
		# Check whether to add /1 or /2 to the read identifier.
		readName = "@"+bamRead.query_name+"/2"
		if(bamRead.is_read1):
			readName = "@"+bamRead.query_name+"/1"
		
		symbolQualities = ''.join([chr((x+33)) for x in bamRead.get_forward_qualities()])	# Convert the Q-Score qualities to ASCII symbols as in the original FastQ file.
		fqEntry = readName + "\n" + bamRead.get_forward_sequence() + "\n+\n" + symbolQualities
		return fqEntry
	
	
	
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
				varcoFile.write("Variant\tSample\tChrom\tStart\tStop\n")
				
				# Iterate over all the variants and their contexts.
				for variant, varContext in variantContextMap.items():
					varcoFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+  varContext[0] +"\t"+ str(varContext[1]) +"\t"+ str(varContext[2]) +"\n")
			self.vaseLogger.info("Finished writing variants and their contexts to " +varConOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write to " +varConOutPath)
			exit()
	
	
	# Writes the VCF variants and their contexts to a separate file.
	def writeVariantsContextsGzip(self, variantContextMap, variantSampleMap, varConOutPath):
		try:
			# Write the variants with their contexts to a specified output file.
			self.vaseLogger.info("Start writing variants and their contexts to " +varConOutPath)
			with gzip.open(varConOutPath, 'wt') as varcoFile:
				varcoFile.write("Variant\tSample\tChrom\tStart\tStop\n")
				
				# Iterate over all the variants and their contexts.
				for variant, varContext in variantContextMap.items():
					varcoFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+ varContext[0] +"\t"+ str(varContext[1]) +"\t"+ str(varContext[2]) +"\n")
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
					try:
						bamFile = pysam.AlignmentFile(variantBamFileMap[variant], 'rb')
						varBreadFile.write(variant +"\t"+ variantSampleMap[variant] +"\t"+ "\t".join([str(x.query_name) for x in bamReads]) +"\n")
						bamFile.close()
					except IOError as ioe:
						self.vaseLogger.warning("Could not open BAM file " +variantBamFileMap[variant]+ " to obtain read information")
			
			self.vaseLogger.info("Finished writing variants and their associated BAM reads to " +varBreadOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write to " +varBreadOutPath)
			exit()
	
	
	
	#Writes the NIST BAM reads associated with variants to a specified output file.
	def writeNistVariantBamReads(self, nistVariantReadMap, nistBreadOutPath):
		try:
			with open(nistBreadOutPath, 'w') as nistBreadFile:
				nistBreadFile.write("Variant\tReads\n")
				
				for variant, nistReads in nistVariantReadMap.items():
					nistBreadFile.write(variant + "\t" + "\t".join([str(x.query_name) for x in nistReads]) + "\n")
				nistBreadFile.close()
		except IOError as ioe:
			self.vaseLogger.critical("Could not write NIST variants read to " +ioe.filename)
			exit()
	
	
	#Writes the NIST BAM reads associated with variants to a specified gzip output file.
	def writeNistVariantBamReadsGzip(self, nistVariantReadMap, nistBreadOutPath):
		try:
			with gzip.open(nistBreadOutPath, 'wt') as nistBreadFile:
				nistBreadFile.write("Variant\tReads\n")
				
				for variant, nistReads in nistVariantReadMap.items():
					nistBreadFile.write(variant + "\t" + "\t".join([str(x.query_name) for x in nistReads]) + "\n")
				nistBreadFile.close()
		except IOError as ioe:
			self.vaseLogger.critical("Could not write NIST variants read to " +ioe.filename)
			exit()
	
	
	
	# Returns the name for the fastq out file.
	def setFastqOutPath(self, outPath, fR):
		valName = outPath+ "/nistVaSe_" +str(datetime.now().date())+ "_R2.fastq.gz"
		if(fR=="F"):
			valName = outPath+ "/nistVaSe_" +str(datetime.now().date())+ "_R1.fastq.gz"
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
