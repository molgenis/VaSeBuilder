import logging
from VaSeUtilHelper import VaSeUtilHelper

class LogInfo:
	def __init__(self, vaseuhelper):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = vaseuhelper
	
	
	# Runs the loginfo program and prints info from the log file.
	def main(self, vaseLogLoc, logFilter):
		self.vaseUtilLogger.info("Running VaSe util LogInfo")
		self.processLogFile(vaseLogLoc, logFilter)
		self.vaseUtilLogger.info("Finished running VaSe util LogInfo")
	
	
	# Processes the log file and prints all lines or lines satisfying the filter
	def processLogFile(self, vaseLogLoc, logFilter):
		try:
			with open(vaseLogLoc, 'r') as vaseLogFile:
				for fileLine in vaseLogFile:
					fileLineElements = fileLine.split(' ')
					if(self.vuh.passes_filter(fileLineElements[3], logFilter)):
						print(fileLine.strip())
		except IOError as ioe:
			self.vaseUtilLogger.warning("Could not open log file " +str(ioe.filename))
