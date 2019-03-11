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
	
	
	# Checks whether a provided output location is valid (does the folder exist, if the output path is a file then use the folder name, etc)
	def checkOutputLocation(self, outLoc):
		if(os.path.isfile(outLoc)):
			# Get the folder name
		elif(os.path.isdir(outLoc)):
			# Everything is ok
	
	
	# Checks the program parameters
	def checkParameters(self, paramList):
		for param in paramList:
			if(param=="log"):
				
			else:
				# All other parameters read specific files.
				if(not self.checkFileExists(paramList[param])):
					self.vaseEvalLogger.critical("Could not read " +str(param)+ " file " +str(paramList[param]))
					return False
		return True
