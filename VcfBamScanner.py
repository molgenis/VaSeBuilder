import logging

class VcfBamScanner:
	#Constructor that creates two empty hashmaps (dicts)
	def __init__(self):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.vcfSampleMap = {}
		self.bamSampleMap = {}
	
	
	
	#Scans the folders containing VCF files and returns a map that links sample ids with vcf files.
	def scanVcfFolders(self, vcfFolders):
		vcfSampleMap = {}
		self.vaseLogger.info("Start scannig VCF files")
		
		for vcfFolder in vcfFolders:
			self.vaseLogger.info("Scanning VCF files in " +vcfFolder)
			for vcfFileName in os.listdir(vcfFolder):
				if vcfFileName.endswith(".vcf"):
					#Do something to scan the file for the sample
					self.vaseLogger.debug("Scanning VCF file " +vcfFileName)
					vcfFile = pysam.VariantFile(vcfFolder+"\\"+vcfFileName, 'r')
					sampleid = ""
					self.vcfSampleMap[sampleid] = vcfFile
		
		self.vaseLogger.info("Finished scanning VCF files")
		return self.vcfSampleMap
	
	
	
	#Scans the folders containing BAM files and returns a map that links sample ids with bam files.
	def scanBamFolders(self, bamFolders):
		bamSampleMap = {}
		self.vaseLogger.info("Start scanning BAM files")
		
		for bamFolder in bamFolders:
			self.vaseLogger.info("Scanning BAM files in " +bamFolder)
			for bamFileName in os.listdir(bamFolder):
				if(bamFileName.endswith(".bam")):
					self.vaseLogger.debug("Scanning BAM file " +bamFileName)
					bamFile = pysam.AlignmentFile(bamFolder+"\\"+bamFileName, 'rb')
					sampleid = bamFile.header["RG"][0]["SM"]	#The sample identifier
					self.bamSampleMap[sampleid] = bamFile
		
		self.vaseLogger.info("Finished scanning BAM files")
		return self.bamSampleMap
	
	
	
	#Returns the map that links VCF files and samples
	def getVcfSampleMap(self):
		return self.vcfSampleMap
	
	
	
	#Returns the map that links BAM files and samples
	def getBamSampleMap(self):
		return self.bamSampleMap
	
	
	
	#Returns a map that links VCF files to BAM files
	def getVcfToBamMap(self):
		vcfToBamMap = {}
		for sampleid in self.vcfSampleMap:
			vcfToBamMap[self.vcfSampleMap[sampleid]] = self.bamSampleMap[vcfSampleMap]
		return vcfToBamMap
