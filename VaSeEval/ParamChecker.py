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
			self.vaseLogger.debug("File " +fileLoc+ " exists")
			return True
		self.vaseLogger.debug("File " +fileLoc+ " does not exist")
		return False
	
	
	# Returns whether a folder exists or not
	def folderExists(self, folderLoc):
		return os.path.isdir(folderLoc)
	
	
	# Returns a proper output location to use to write the output files to.
	def getProperOutLocation(self, outLoc):
		if(os.path.isfile(outLoc)):
			if(self.folderExists(os.path.dirname(outLoc))):
				return os.path.dirname(outLoc)
			else
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
			if(param=="log"):
				if(not self.isValidOutputLocation(paramList[param])):
					self.vaseEvalLogger.warning("Provided log location is invalid. Will use default output location")
			elif(param=="out"):
				if(not (self.isValidOutputLocation(paramList[param]))):
					return False
			else:
				# All other parameters read specific files.
				if(not self.checkFileExists(paramList[param])):
					self.vaseEvalLogger.critical("Could not read " +str(param)+ " file " +str(paramList[param]))
					return False
		return True
	
	
	def getParameterMap(self, paramList):
