import logging
from datetime import datetime
from Bio import SeqIO
import gzip
import pysam
import numpy


class VaSeBuilder:
	#Constructor that saves the identifier, date and time of the current 
	def __init__(self, vaseId):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.creationId = vaseId
		self.creationDate = datetime.now().date()
		self.creationTime = datetime.now().time()
		self.variantContextMap = {}
		self.variantBamReadMap = {}
		self.nistVariantReadMap = {}
		self.vaseLogger.info("VaSeBuilder: " +self.creationId+ " ; " +str(self.creationDate)+ " ; " +str(self.creationTime))
	
	
	
	#Creates the new FastQ validation dataset by replacing NIST reads containing a VCF variant with BAM reads from patients. Returns true at the end to indicate the process is done.
	def buildValidationSet(self, vcfSampleMap, bamSampleMap, nistBamLoc, vcfBamMap, fastqFPath, fastqRPath, fastqOutPath, varConOutPath, varBreadOutPath, nistBreadOutPath):
		try:
			nistBam = pysam.AlignmentFile(nistBamLoc, 'rb')
			self.vaseLogger.info("Start building the valdation set")
			
			#Iterate over the samples to use for building the validation set.
			for sampleId in vcfSampleMap:
				#Check in the VCF BAM link map that there is a BAM file for the sample as well.
				if(sampleId in vcfBamMap):
					try:
						vcfFile = pysam.VariantFile(vcfSampleMap[sampleId], 'r')
						
						#Loop over the variants in the VCF file. Prior to identifying the BAM reads, it is first checked whether the variant is in a previously established 
						for vcfVar in vcfFile.fetch():
							variantId = self.getVcfVariantId(vcfVar)
							
							#Get the BAM reads fr the variant and determine the variant context.
							if(not self.isInContext(vcfVar.chrom, vcfVar.pos)):
								self.variantBamReadMap[variantId] = self.getVariantReads(vcfVar.chrom, vcfVar.pos, bamSampleMap[sampleId])	#Obtain all patient BAM reads containing the VCF variant and their read mate.
								self.variantContextMap[variantId] = self.determineContext(self.variantBamReadMap[variantId])	#Save the context start and stop for the variant in the VCF variant context map.
								self.nistVariantReadMap[variantId] = self.getVariantReads(vcfVar.pos, nistBam)	#Obtain all NIST BAM reads containing the VCF variant. and their read mate.
						vcfFile.close()
					except IOError as ioe:
						self.vaseLogger.warning("Could not establish data for " +sampleId)
			nistBam.close()
			
			#Write data used to build the new FastQ to output files.
			self.writeVariantsContexts(self.variantContextMap, varConOutPath)	#Write the context start and stop for each used variant to a separate file.
			self.writeVariantsBamReads(self.variantBamReadMap, varBreadOutPath)	#Write the associated BAM reads for each used variant to a seperate file.
			self.writeVariantsContexts(self.nistVariantReadMap, nistBreadOutPath)	#Write the associated NIST BAM reads for each used variant to a separate file.
			
			#Obtain a list of NIST reads to skip when iterating over the NIST FastQ.
			nistReadsToSkip = numpy.array(list(self.nistVariantReadMap.values())).ravel()	#Set up a list of all NIST reads to skip.
			
			#Make the new FastQ files that can be used to run in the NGS_DNA pipeline along real sample data
			self.buildFastQ(fastqFPath, self.setFastqOutPath(fastqOutPath), nistReadsToSkip, self.variantBamReadMap, 'F')	#Build the R1 fastq file.
			self.buildFastQ(fastqRPath, self.setFastqOutPath(fastqOutPath), nistReadsToSkip, self.variantBamReadMap, 'R')	#Build the R2 fastq file.
			self.vaseLogger.info("Finished building the validation set")
		except:
			self.vaseLogger.critical("Could not open the NIST BAM and can therefore not build the new FastQ files.")
			exit()
	
	
	
	#Returns the BAM reads containing the specific vcf variant as well as their read mate.
	def getVariantReads(self, vcfVariantChr, vcfVariantPos, bamFileLoc):
		try:
			bamFile = pysam.AlignmentFile(bamFileLoc, 'rb')
			
			#Obtain all the variant reads overlapping with the variant and their mate reads.
			variantReads = []
			for vread in bamFile.fetch(vcfVariantChr, vcfVariantPos-1, vcfVariantPos+1):
				variantReads.append(vread)	#Add the BAM read to the list of reads.
				variantReads.append(bamFile.mate(vread))	#Add the mate of the current BAM read to the list as well
			bamFile.close()
			
			return variantReads
		except IOError as ioe:
			self.vaseLogger.warning("Could not obtain BAM reads from " +bamFileLoc)
			return None
	
	
	
	#Determines the start and stops of the variant context (please see the documentation for more information).
	def determineContext(self, bamVariantReads):
		#First determine the context start by sorting the reads on leftmost position in ascending order.
		bamVariantReads.sort(key=lambda x:x.reference_start, reverse=False)
		contextChrom = bamVariantsReads[0].reference_name
		contextStart = bamVariantReads[0].reference_start
		
		#Second determine the context stop by iterating over the reads and calculating the rightmost position of the reads.
		contextStop = 0
		for bvRead in bamVariantReads:
			stopPos = bvRead.reference_start + bvRead.infer_read_length()
			if(stopPos > contextStop):
				contextStop = stopPos
		
		return [contextChrom, contextStart, contextStop]
	
	
	
	#Returns whether a certain variant is in an already established variant context.
	def isInContext(self, vcfVarChrom, vcfVarPos):
		for vcfVar, context in self.variantContextMap:
			if(vcfVarChrom == context[0]):
				if(vcfVarPos >= context[1] and vcfVarPos <= context[2]):
					return True
		return False
	
	
	
	#Returns whether a variant has already been used before in the analysis.
	def variantAlreadyProcessed(self, vcfVarId):
		return vcfVarId in self.variantContextMap
	
	
	
	#Builds the new FastQ file to be used.
	def buildFastQ(self, nistFastqIn, fastqOutPath, nistReadsToSkip, patientBamReads, fR):
		#Build the new NIST FastQ file replacing some of its own reads with patient reads containing VCF variants.
		try:
			#Write NIST fastq entries to the new FastQ file but exclude NIST reads overlapping with the position of a VCF variant.
			with gzip.open(nistFastqIn, "rt") as nistFq, gzip.open(fastqOutPath, 'wt') as newNistFq:
				for fqRead in SeqIO.parse(nistFq, "fastq"):
					if(fqRead.id not in nistFqReads):
						newNistFq.write(fqRead.format("fastq"))
				
				#Add the patient BAM reads containing a VCF variant to the new FastQ file.
				for vcfvar in patientBamReads:
					for bamRead in patientBamReads:
						#Check if the BAM read is R1 or R2.
						if(self.isRequiredRead(bamRead, fR)):
							newNistFq.write(getBamReadAsFastQ(bamRead)+"\n")
		
		#Throw an exception if one of the files does not exist.
		except IOError as ioe:
			if(ioe.filename==nistFastqIn):
				self.vaseLogger.critical("The supplied NIST FastQ file could not be found.")
			if(ioe.filename==fastqOutPath):
				self.vaseLogger.critical("A FastQ file could not be written to the provided output location.")
			exit()
	
	
	
	#Checks if a read is read 1 (R1) or read 2 (R2).
	def isRequiredRead(self, bamRead, fR):
		if(fR=="F"):
			if(bamRead.is_read1):
				return True
		else:
			if(bamRead.is_read2):
				return True
		return False
	
	
	
	#Obtains all required info from the 
	def getBamReadAsFastQ(self, bamRead):
		fqEntry = bamRead.id + "\n" + bamRead.get_forward_sequence() + "\n+" + bamRead.get_forward_qualities()
	
	
	
	#Returns the identifier of the current VaSeBuilder object.
	def getCreationId(self):
		return self.creationId
	
	
	
	#Returns the date the current VaSeBuilder object has been made.
	def getCreationDate(self):
		return self.creationDate
	
	
	
	#Returns the time the current VaSeBuilder object has been made.
	def getCreationTime(self):
		return self.creationTime
	
	
	
	#Returns all variant contexts.
	def getVariantContexts(self):
		return self.variantContextMap
	
	
	
	#Returns the context start and stop for a specified VCF variant.
	def getVariantContext(self, vcfVariant):
		if(vcfVariant in self.variantContextMap):
			return self.variantContextMap[vcfVariant]
		else:
			return None
	
	
	
	#Adds a VCF variant and its context start and stop if it isn't hasn't been established before.
	def addVariantContext(self, vcfVariant, contextStart, contextStop):
		if(vcfVariant not in self.variantContextMap):
			self.variantContextMap[vcfVariant] = [contextStart, contextStop]
	
	
	
	#Writes the VCF variants and their contexts to a separate file.
	def writeVariantsContexts(self, variantContextMap, varConOutPath):
		try:
			#Write the variants with their contexts to a specified output file.
			self.vaseLogger.info("Start writing variants and their contexts to " +varConOutPath)
			with open(varConOutPath) as varcoFile:
				varcoFile.write("Variant\tChrom\tStart\tStop\n")
				
				#Iterate over all the variants and their contexts.
				for variant, varContext in variantContextMap:
					varcoFile.write(variant +"\t"+ varContext[0] +"\t"+ varContext[1] +"\t"+ varContext[2] +"\n")
			self.vaseLogger.info("Finished writing variants and their contexts to " +varConOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write to " +varConOutPath)
			exit()
	
	
	
	#Writes the VCF variants and their associated BAM reads to a separate file.
	def writeVariantsBamReads(self, variantBamReadMap, varBreadOutPath):
		try:
			#Write the variants with their contexts to a specified output file.
			self.vaseLogger.info("Start writing variants and their associated BAM reads to " +varBreadOutPath)
			with open(varBreadOutPath) as varBreadFile:
				varBreadFile.write("Variant\tStart\tStop\n")
				
				#Iterate over all the variants and their contexts.
				for variant, bamReads in variantBamReadMap:
					varBreadFile.write(variant +"\t"+ "\t".join([str(x.id) for x in bamReads]) +"\n")
			self.vaseLogger.info("Finished writing variants and their associated BAM reads to " +varBreadOutPath)
		
		except IOError as ioe:
			self.vaseLogger.critical("Could not write to " +varBreadOutPath)
			exit()
	
	
	
	#Returns the name for the fastq out file.
	def setFastqOutPath(self, outPath, fR):
		valName = outPath+ "/nistVaSe_" +datetime.now().date()+ "_" +datetime.now().time()+ "_R2.fastq.gz"
		if(fR=="F"):
			valName = outPath+ "/nistVaSe_" +datetime.now().date()+ "_" +datetime.now().time()+ "_R1.fastq.gz"
		return valName
	
	
	
	#Returns an identifier for a VCF variant. If the identifier is '.' then one will be constructed as 'chrom_pos'. 
	def getVcfVariantId(self, vcfVariant):
		variantId = vcfVariant.id
		if(variantId=='.' or variantId==None):
			variantId = vcfVariant.chrom +"_"+ vcfVariant.pos
		return variantId
