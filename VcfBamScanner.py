#!/usr/bin/env python
import logging
import os
import pysam

class VcfBamScanner:
	# Constructor that creates two empty hashmaps (dicts)
	def __init__(self):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.vcfSampleMap = {}
		self.bamSampleMap = {}
	
	
	
	# Scans the folders containing VCF files and returns a map that links sample ids with vcf files.
	def scanVcfFolders(self, vcfFolders):
		self.vaseLogger.info("Start scannig VCF files")
		
		for vcfFolder in vcfFolders:
			if(os.path.isdir(vcfFolder)):
				self.vaseLogger.info("Scanning VCF files in " +vcfFolder)
				for vcfFileName in os.listdir(vcfFolder):
					if(vcfFileName.endswith(".vcf") or vcfFileName.endswith(".vcf.gz")):
						# Do something to scan the file for the sample
						try:
							vcfFile = pysam.VariantFile(vcfFolder+"/"+vcfFileName, 'r')
							self.vaseLogger.debug("Scanning VCF file " +vcfFolder+ "/" +vcfFileName)
							if(len(vcfFile.header.samples) > 0):
								sampleid = vcfFile.header.samples[0]	# This is the sample identifier
								self.vcfSampleMap[sampleid] = vcfFolder+"/"+vcfFileName
							vcfFile.close()
						except IOError as ioe:
							self.vaseLogger.warning("VCF file " +vcfFolder+ "/" +vcfFileName+ " does not exist")
			else:
				self.vaseLogger.info("Folder " +vcfFolder+ " does not exist.")
		self.vaseLogger.info("Finished scanning VCF files")
		return self.vcfSampleMap
	
	
	
	# Scans the folders containing BAM files and returns a map that links sample ids with bam files.
	def scanBamFolders(self, bamFolders):
		self.vaseLogger.info("Start scanning BAM files")
		
		# Scan BAM files in all provided folders
		for bamFolder in bamFolders:
			if(os.path.isdir(bamFolder)):
				self.vaseLogger.info("Scanning BAM files in " +bamFolder)
				for bamFileName in os.listdir(bamFolder):
					if(bamFileName.endswith(".bam")):
						try:
							bamFile = pysam.AlignmentFile(bamFolder+"/"+bamFileName, 'rb')
							self.vaseLogger.debug("Scanning BAM file " +bamFolder+ "/" +bamFileName)
							if(self.bamHasSampleName(bamFile)):
								sampleid = bamFile.header["RG"][0]["SM"]	# The sample identifier
								self.bamSampleMap[sampleid] = bamFolder+"/"+bamFileName
							bamFile.close()
						except IOError as ioe:
							self.vaseLogger.warning("BAM file " +bamFolder+ "/" +bamFileName+ " does not exist")
			else:
				self.vaseLogger.info("Folder " +bamFolder+ " does not exist.")
		self.vaseLogger.info("Finished scanning BAM files")
		return self.bamSampleMap
	
	
	# Checks whether the BAM file contains a sample name.
	def bamHasSampleName(self, bamFile):
		if("RG" in bamFile.header):
			if(len(bamFile.header["RG"]) > 0):
				if("SM" in bamFile.header["RG"][0]):
					return True
		return False
	
	
	
	# Returns the map that links VCF files and samples
	def getVcfSampleMap(self):
		return self.vcfSampleMap
	
	
	
	# Returns the map that links BAM files and samples
	def getBamSampleMap(self):
		return self.bamSampleMap
	
	
	
	# Returns a map that links VCF files to BAM files
	def getVcfToBamMap(self):
		vcfToBamMap = {}
		for sampleid in self.vcfSampleMap:
			if(sampleid in self.bamSampleMap):
				vcfToBamMap[self.vcfSampleMap[sampleid]] = self.bamSampleMap[sampleid]
		return vcfToBamMap
