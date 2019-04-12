#!/usr/bin/env python
import logging
import os

class ParamChecker:
	# Constructor
	def __init__(self):
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
	
	
	# Checks whether a file exists.
	def checkFileExists(self, fileLoc):
		if(os.path.isfile(fileLoc)):
			self.vaseEvalLogger.debug("File " +fileLoc+ " exists")
			return True
		self.vaseEvalLogger.debug("File " +fileLoc+ " does not exist")
		return False
	
	
	# Returns whether a folder exists or not
	def folderExists(self, folderLoc):
		return os.path.isdir(folderLoc)
	
	
	# Returns a proper output location to use to write the output files to.
	def getProperOutLocation(self, outLoc):
		if(os.path.isfile(outLoc)):
			if(self.folderExists(os.path.dirname(outLoc))):
				return os.path.dirname(outLoc)
			else:
				return ""
		elif(os.path.isdir(outLoc)):
			return outLoc
		else:
			return ""
	
	
	# Checks whether a provided output location is valid (does the folder exist, if the output path is a file then use the folder name, etc)
	def isValidOutputLocation(self, outLoc):
		if(not (os.path.isfile(outLoc))):
			return self.folderExists(outLoc)
		else:
			return False
	
	
	# Checks the program parameters
	def checkParameters(self, paramList):
		for param in paramList:
			if(param=="out"):
				if(not (self.isValidOutputLocation(paramList[param]))):
					self.vaseEvalLogger.info("Output location is invalid")
					return False
			else:
				if(param != 'log'):
					# All other parameters read specific files.
					if(not self.checkFileExists(paramList[param])):
						self.vaseEvalLogger.critical("Could not read " +str(param)+ " file " +str(paramList[param]))
						return False
		return True
	
	
	
	# Check the logging parameter to determine where to write the logfile to.
	def checkLog(self, logParam):
		# Check the filepath for the log output file
		if(logParam == None):
			self.logLocation = "VaSeEvaluator.log"
			return "VaSeEvaluator.log"
		
		else:
			logloc = ""
			
			# Check the location of the log file if the --log parameter has been set.
			if( not(os.path.isfile(logParam)) and (logParam.endswith(".log") or logParam.endswith(".txt")) ):
				logloc = logParam
			else:
				logloc = "VaSeEvaluator.log"
			
			# Check to make sure the provided --log parameter value is not a directory (Directories could be named 'something.log')
			if(os.path.isdir(logParam)):
				logloc = logParam + "/VaSeEvaluator.log"
			self.logLocation = logloc
			return self.logLocation
