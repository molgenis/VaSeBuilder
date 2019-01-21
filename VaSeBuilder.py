import logging
from datetime import datetime


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
	def buildValidationSet(self, vcfSampleMap, bamSampleMap, nistBamLoc, fastqOutPath):
		nistBam = pysam.AlignmentFile(nistBamLoc, 'rb')
		self.vaseLogger.info("Start building the valdation set")
		
		for sampleId in vcfSampleMap:
			#Loop over the variants in the VCF file. Prior to identifying the BAM reads, it is first checked whether the variant is in a previously established 
			for vcfVar in vcfSampleMap[sampleId].fetch():
				if(not self.isInContext(vcfVar.id)):
					self.variantBamReadMap[vcfVar.id] = self.getVariantReads(vcfVar.pos, bamSampleMap[sampleId])	#Obtain all patient BAM reads containing the VCF variant and their read mate.
					self.variantContextMap[vcfVar.id] = self.determineContext(vcfVarBamReads)	#Save the context start and stop for the variant in the VCF variant context map.
					self.nistVariantReadMap[vcfVar.id] = self.getVariantReads(vcfVarPos, nistBam)	#Obtain all NIST BAM reads containing the VCF variant. and their read mate.
		
		#Write the context start and stop for each used variant to a separate file (so it can be viewed afterwards).
		self.writeVariantsContexts()
		
		#Make the new FastQ file that can be used to run in the NGS_DNA pipeline along real sample data
		self.buildFastQ(fastqOutPath)
		self.vaseLogger.info("Finished building the validation set")
		return True
	
	
	
	#Returns the BAM reads containing the specific vcf variant as well as their read mate.
	def getVariantReads(self, vcfVariantPos, bamFile):
		variantReads = bamFile.fetch(vcfVariantPos-1, vcfVariantPos+1)
		
		#Obtain the read mate for each variant read. It will first be saved in a separate list just to be safe.
		variantMateReads = []
		for variantRead in variantReads:
			variantMateReads.append(variantRead.mate())
		
		#Check if the number of mate reads equals the number of variant reads.
		if(len(variantMateReads) < len(variantReads)):
			self.vaseLogger.debug("Fewer mate reads than variant reads")
		elif(len(varianMatetReads) > len(variantReads)):
			self.vaseLogger.debug("More mate reads than variant reads")
		else:
			self.vaseLogger.debug("Number of read mates equals number of variant reads")
		
		#Combine the list of variant reads and 
		allReads = variantReads + variantMateReads
		return allReads
	
	
	
	#Determines the start and stops of the variant context (please see the documentation for more information).
	def determineContext(self, bamVariantReads):
		#First determine the context start by sorting the reads on leftmost position in ascending order.
		bamVariantReads.sort(key=lambda x:x.reference_start, reverse=False)
		contextStart = bamVariantReads[0].reference_start
		
		#Second determine the context stop by iterating over the reads and calculating the rightmost position of the reads.
		contextStop = 0
		for bvRead in bamVariantReads:
			stopPos = bvRead.reference_start + bvRead.infer_read_length()
			if(stopPos > contextStop):
				contextStop = stopPos
		
		return [contextStart, contextStop]
	
	
	
	#Checks whether a VCF variant is located within an earlier established variant context.
	def isInContext(self, vcfVariant):
		return (vcfVariant in self.variantContextMap)
	
	
	
	#Builds the new FastQ file to be used.
	def buildFastQ(self, outFilePath):
		print("Implementing...")
	
	
	
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
	def writeVariantContexts(self, varConOutPath):
		variantContextOutFile = open(varConOutPath, 'w')
		self.vaseLogger.info("Start writing VCF variant contexts to " + varConOutPath)
		variantContextOutFile.write("ID\tStart\tStop\n")
		
		for vcfVariant in self.variantContextMap:
			if(vcfVariant in self.variantContextMap):
				variantContextOutFile.write(vcfVariant + "\t" + self.variantContextMap[vcfVariant][0] + "\t" + self.variantContextMap[vcfVariant][1] + "\n")
		variantContextOutFile.close()
		self.vaseLogger.info("Finished writing VCF variant context data")
