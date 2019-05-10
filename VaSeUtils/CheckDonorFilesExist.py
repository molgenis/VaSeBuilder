import os
import logging

class CheckDonorFilesExist:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Checks whether the donor VCF/BAM files still exist at the given location.
	def main(donorListFile, samplesToCheck=None):
		self.vaseUtilLogger.info("Run VaSe util CheckDonorFilesExist")
		try:
			with open(donorListFile, 'r') as dlFile:
				next(dlFile)	# Skip the header line of the file
				for fileLine in dlFile:
					fileLineData = fileLine.strip().split("\t")
					donorFiles = fileLineData[1].split(" ; ")
					
					# Check whether the used donor files exist.
					if(samplesToCheck is not None):
						if(fileLineData[0] in samplesToCheck):
							for dfile in donorFiles:
								if(not os.path.isfile(dfile)):
									self.vaseUtilLogger.info("Donor file " +str(dfile)+ " could not be found. Maybe it has been moved, renamed or deleted?")
		except IOError as ioe:
			self.vaseUtilLogger.info("File " +str(ioe.filename)+ " containing used donor files does not seem to exist.")
			exit()
		self.vaseUtilLogger.info("Finished running VaSe util CheckDonorFilesExist")
