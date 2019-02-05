#!/usr/bin/env python
import logging
import os

class ParamChecker:
	# Constructor that creates two empty arrays that 
	def __init__(self):
		self.vaseLogger = logging.getLogger("VaSe_Logger")
		self.vcfFolders = []
		self.bamFolders = []
		self.nistBam= ""
		self.fastqIn1 = ""
		self.fastqIn2 = ""
		self.fastqOutLocation = ""
		self.varConOutLocation = ""
		self.varBreadOutLocation = ""
		self.nistBreadOutLocation = ""
		self.logLocation = ""
	
	
	
	# Check the logging parameter to determine where to write the logfile to.
	def checkLog(self, logParam):
		# Check the filepath for the log output file
		if(logParam == None):
			self.logLocation = "VaSeBuilder.log"
			return "VaSeBuilder.log"
		
		else:
			logloc = ""
			
			# Check the location of the log file if the --log parameter has been set.
			if( not(os.path.isfile()) and (logParam.endswith(".log") or logParam.endswith(".txt")) ):
				logloc = logParam
			else:
				logloc = "VaSeBuilder.log"
			
			# Check to make sure the provided --log parameter value is not a directory (Directories could be named 'something.log')
			if(os.path.isdir(logParam)):
				logloc = logParam + "/VaSeBuilder.log"
			self.logLocation = logloc
			return self.logLocation
	
	
	
	# Checks whether provided folders exist
	def checkFoldersExist(self, paramVals, fileExt):
		existingFolders = []
		
		# Check whether the provided folders exist.
		for parval in paramVals:
			foldername = self.getFolderName(parval)
			if(foldername != parval):
				self.vaseLogger.warning("File instead of folder was supplied. Using the folder " +foldername+ " as input folder")
			
			# Check if the supplied value is a folder or not and contains any vcf/bam files.
			if(not (os.path.isdir(foldername))):
				self.vaseLogger.warning("Folder " +foldername+ " was not found and will therefore be skipped")
			else:
				if(self.checkFolderContents(foldername, fileExt)==0):
					self.vaseLogger.warning("Folder " +foldername+ " exists but contains no " +fileExt+ " files")
				else:
					self.vaseLogger.info("Folder " +foldername+ " will be included")
					existingFolders.append(foldername)
		return existingFolders
	
	
	
	# Checks whether at least one file with a provided extension (.vcf or .bam) is present.
	def checkFolderContents(self, folderToCheck, fileExt):
		vbCnt = 0
		for vbFile in os.listdir(folderToCheck):
			if(vbFile.endswith("."+fileExt)):
				vbCnt += 1
		self.vaseLogger.debug("Folder " +folderToCheck+ " contains " +str(vbCnt)+ " " +fileExt+ " files")
		return vbCnt
	
	
	
	# Checks whether a provided file exists.
	def checkFileExist(self, fileLoc):
		if(os.path.isfile(fileLoc)):
			self.vaseLogger.debug("File " +fileLoc+ " exists")
			return True
		self.vaseLogger.debug("File " +fileLoc+ " does not exist")
		return False
	
	
	# 
	def isValidOutputLocation(self, outFileName):
		return (os.path.isdir(os.path.dirname(outFileName)))
	
	
	
	# Checks whether the values of the parameters are correct (do files/folders exist for example)
	# [Function should perhaps be split into smaller functions]
	def checkParameters(self, vaseArgVals):
		
		# Loop over the provided parameters.
		for param in vaseArgVals:
			
			# If the current parameter is vcfin, check whether there are any valid VCF folders to use.
			if(param=="vcfin"):
				vcfFolders = self.checkFoldersExist(vaseArgVals[param], "vcf")
				if(len(vcfFolders)==0):
					self.vaseLogger.critical("No folders containing VCF files were found. Please supply existing folders next time :)")
					return False
				self.vcfFolders = vcfFolders
			
			# If the current parameter is bamin, check whether there are any valid BAM folders to use.
			if(param=="bamin"):
				bamFolders = self.checkFoldersExist(vaseArgVals[param], "bam")
				if(len(bamFolders)==0):
					self.vaseLogger.critical("No folders containing BAM files were found. Please supply existing folders next time :)")
					return False
				self.bamFolders = bamFolders
			
			# If the current parameter is bam, check whether a valid BAM file is provided.
			if(param=="valbam"):
				if(not self.checkFileExist(vaseArgVals[param])):
					self.vaseLogger.critical("No valid NIST BAM file supplied :(")
					return False
				self.nistBam = vaseArgVals[param]
			
			# If the current parameter is valfastq1, check whether a valid R1 fastq file is provided.
			if(param=="valfastq1"):
				if(not self.checkFileExist(vaseArgVals[param])):
					self.vaseLogger.critical("Provided R1 FastQ input file does not exist")
					return False
				self.fastqIn1 = vaseArgVals[param]
			
			# If the current parameter is valfastq2, check whether a valid R2 fastq file is provided.
			if(param=="valfastq2"):
				if(not self.checkFileExist(vaseArgVals[param])):
					self.vaseLogger.critical("Provided R2 FastQ input file does not exist")
					return False
				self.fastqIn2 = vaseArgVals[param]
			
			# If the current parameters is fastqout
			if(param=="fastqout"):
				if(not (os.path.isdir(vaseArgVals[param]))):
					return False
				self.fastqOutLocation = vaseArgVals[param]
			
			# If the current parameter is varcon, check whether a valid output location is provided
			if(param=="varcon"):
				if(not(self.isValidOutputLocation(vaseArgVals[param]))):
					return False
				self.varConOutLocation = vaseArgVals[param]
			
			# If the current parameters is varbread, check whether a valid output location is provided
			if(param=="varbread"):
				if(not(self.isValidOutputLocation(vaseArgVals[param]))):
					return False
				self.varBreadOutLocation = vaseArgVals[param]
			
			# If the current parameter is nistbread, check whether a valid output location is provided
			if(param=="nistbread"):
				if(not(self.isValidOutputLocation(vaseArgVals[param]))):
					return False
				self.nistBreadOutLocation = vaseArgVals[param]
			
		# Return the lists of valid VCF and BAM folders that can be used by the program.
		return True
	
	
	
	# Returns thename of the folder name of a parameter value (if the parameter value is )
	def getFolderName(self, foldername):
		if(os.path.isfile(foldername)):
			return os.path.dirname(foldername)
		return foldername
	
	
	
	# Returns the list of valid VCF folders.
	def getValidVcfFolders(self):
		return self.vcfFolders
	
	
	
	# Returns the list of valid BAM folders.
	def getValidBamFolders(self):
		return self.bamFolders
	
	
	
	# Returns the location of the  NIST BAM file.
	def getNistBam(self):
		return self.nistBam
	
	
	
	# Returns the location and name of the first (R1) fastq input file.
	def getFirstFastqInLocation(self):
		return self.fastqIn1
	
	
	
	# Returns the location and name of the second (R2) fastq input file.
	def getSecondFastqInLocation(self):
		return self.fastqIn2
	
	
	
	# Returns the location(s) and names of the two (R1 and R2) fastq input files.
	def getFastqInLocations(self):
		return [self.fastqIn1, self.fastqIn2]
	
	
	
	# Returns the location of the FastQ file that will be produced by VaSeBuilder.
	def getFastqOutLocation(self):
		return self.fastqOutLocation
	
	
	
	# Returns the location of file that will contain the variants and their context start and stops.
	def getVariantContextOutLocation(self):
		return self.varConOutLocation
	
	
	
	# Returns the location of the file that will contain the variants and their associatied patient BAM reads.
	def getVariantBamReadOutLocation(self):
		return self.varBreadOutLocation
	
	
	
	# Returns the location of the file that will containt the variants and their associated NIST reads.
	def getNistBamReadOutLocation(self):
		return self.nistBreadOutLocation
	
	
	
	# Retuns the location to write the log file(s) to.
	def getLogFileLocation(self):
		return self.logLocation
