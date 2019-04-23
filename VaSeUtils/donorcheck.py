import logging
import gzip

class DonorCheck:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaseUtil_Logger")
	
	
	# Performs the main analysis
	def main(self, readsList, vaseFq1, vaseFq2):
		self.vaseUtilLogger.info("Running VaSe util DonorCheck")
		
		self.vaseUtilLogger.info("Checking the R1 VaSe FastQ file.")
		r1Added = self.checkDonorReadsAdded(vaseFq1, readsList)	# Check if the donor reads have indeed been added to the VaSe R1 FastQ
		self.vaseUtilLogger.info("Added " +str(r1Added)+ " of " str(len(readsList))+ " to the R1 VaSe FastQ file.")
		
		self.vaseUtilLogger.info("Checking the R2 VaSe FastQ file.")
		r2Added = self.checkDonorReadsAdded(vaseFq2, readsList)	# Check if the donor reads have indeed been added to the VaSe R2 FastQ
		self.vaseUtilLogger.info("Added " +str(r2Added)+ " of " +str(len(readsList))+ " to the R2 VaSe FastQ file.")
		
		if(r1Added==len(readsList) and r2Added==len(readsList)):
			self.vaseUtilLogger.info("All donor reads have been added to the VaSe FastQ files")
		else:
			self.vaseUtilLogger.info("Not all donor reads have been added to the VaSe FastQ files")
		self.vaseUtilLogger.info("Finished running VaSe util DonorCheck")
	
	
	# Returns the number of added reads to R1/R2 VaSe fastq data.
	def getNumberOfReadsAdded(self, vaseFqFiles, readIdList):
		addedReads = 0
		for vfqFile in vaseFqFiles:
			addedReads = self.checkDonorReadsAdded(vfqFile, readIdList, addedReads)
		return addedReads
	
	
	# Checks whether the list of donor reads are indeed added to a specified fastq file based on read identifier.
	def checkDonorReadsAdded(self, gzResultsFile, donorReadList, addedCount):
		with gzip.open(gzResultsFile, 'rt') as gzFile:
			for fileLine in gzFile:
				fileLine = fileLine.strip()
				if(fileLine.startswith('@')):
					if(fileLine[1:] in donorReadList):
						addedCount = addedCount + 1
					else:
						self.vaseUtilLogger.info("Read " +str(fileLine)+ " was not added.")
					next(gzFile)	# Skip the read sequence line
					next(gzFile)	# Skip the '+' line
					next(gzFile)	# Skip the read qualities line
		return addedCount
