import logging
import gzip

class AcceptorCheck:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vaseUtilLogger.info("Running VaSe util AcceptorCheck")
	
	# Runs the check to determine if none of the variant context acceptor BAM reads are in the two VaSe produced FastQ files.
	def main(self, readsList, vaseFq1, vaseFq2):
		r1Removed = checkAcceptorReadsRemoved(vaseFq1, readsList)	# Check if the acceptor reads have indeed not been placed to the VaSe R1 FastQ
		r2Removed = checkAcceptorReadsRemoved(vaseFq2, readsList)	# Check if the acceptor reads have indeed not been added to R2
		
		self.vaseUtilLogger.info("Excluded " +str(r1Removed)+ " of " +str(len(readsList))+ " acceptor reads from the R1 file.")
		self.vaseUtilLogger.info("Excluded " +str(r2Removed)+ " of " +str(len(readsList))+ " acceptor reads from the R2 file.")
		if(r1Removed==len(readsList) and r2Removed==len(readsList)):
			self.vaseUtilLogger.info("All acceptor reads have been excluded from the new VaSe FastQ files")
		else:
			self.vaseUtilLogger.info("There are still some acceptor reads left in th VaSe FastQ files")	
	
	
	# Checks whether the identified acceptor reads have indeed been removed from a specified VaSe FastQ file.
	def checkAcceptorReadsRemoved(self, gzResultsFile, acceptorReadList):
		removalCount = len(acceptorReadList)
		with gzip.open(gzResultsFile, 'rt') as gzFile:
			for fileLine in gzFile:
				fileLine = fileLine.strip()
				if(fileLine.startswith('@')):
					if(fileLine[1:] in acceptorReadList):
						self.vaseUtilLogger.info("Read " +str(fileLine[1:])+ " was not removed")
						removalCount = removalCount - 1
					next(gzFile)	# Skip the sequence line
					next(gzFile)	# Skip the line with '+'
					next(gxFile)	# Skip the sequence qualities line
		return removalCount
