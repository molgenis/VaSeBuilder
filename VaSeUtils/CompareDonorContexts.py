import logging

class CompareDonorContexts:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	def main(self)
	
	
	# Reads the donor context file
	def readDonorContextFile(self, contextFileLoc):
		try:
			with open(contextFileLoc, 'r') as dContextFile:
				next(dContextFile)	# Skip the header line.
				for fileLine in dContextFile:
					
		except IOError as ioe:
			self.vaseUtilLogger.warning("Could not read donor context file " +str(ioe.filename))