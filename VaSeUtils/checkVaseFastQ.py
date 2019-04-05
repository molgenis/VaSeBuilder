import logging
import gzip
from donorcheck import DonorCheck
from acceptorcheck import AcceptorCheck

class CheckVaSeFastQ:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		
	
	# Checks that the VaSe R1 and R2 FastQ have the correct number of lines.
	def main(self, templatefq1, vasefq1, templatefq2, vasefq2, donorReadList, acceptorReadList):
		self.vaseUtilLogger.info("Running VaSe util CheckVaSeFastQ")
		self.checkVaseFastQ(templatefq1, vasefq1, templatefq2, vasefq2, donorReadList, acceptorReadList)
		self.vaseUtilLogger.inf("Finished running VaSe util CheckVaSeFastQ")
	
	
	# Checks the VaSe produced FastQ files (is the size as should be, are all required template reads added to the VaSe produced fastq files)
	def checkVaseFastQ(tmplFq1, vaseFq1, tmplFq2, vasesFq2, donorReadsFile, acceptorReadsFile):
		donorReads = readDonorReadList(donorReadsFile)
		acceptorReads = readAcceptorReadList(acceptorReadsFile)
		
		# Check the VaSe produced R1 fastq file.
		if(checkVaseFastqSize(tmplFq1, vaseFq1, donorReads, acceptorReads)):
			self.vaseUtilLogger.info("VaSe produced FastQ R1 file " +str(vaseFq1)+ "seems to have the expected size.")
		else:
			self.vaseUtilLogger.info("VaSe produced FastQ R1 file " +str(vaseFq1)+ "does not seem to have the expected size.")
		
		# Check the VaSe produced R2 fastq file.
		if(checkVaseFastqSize(tmplFq2, vaseFq2, donorReads, acceptorReads)):
			self.vaseUtilLogger.info("VaSe produced FastQ R2 file " +str(vaseFq2)+ "seems to have the expected size.")
		else:
			self.vaseUtilLogger.info("VaSe produced FastQ R2 file " +str(vaseFq2)+ "does not seem to have the expected size.")
		
		# Lastly check that the reads of the VaSe FastQ and acceptor reads are the same as those in the template (original) fastq files
		#checkReadsInFastQs()	# Check for R1
		#checkReadsInFastQs()	# Check for R2
	
	
	# Compares an original (template) and VaSe produced fastq file
	def checkVaseFastqSize(tmplFq, vaseFq, donorreads, acceptorreads):
		numOfDonorReads = checkDonorReadsAdded(vaseFq, donorreads)
		numOfAcceptorReads = checkAcceptorReadsRemoved(vaseFq, acceptorreads)
		
		# Get the number of lines for the template/original fastq file
		tmplfqNumOfLines = 0
		if(tmplFq.endswith('.fastq') or tmplFq.endswith('.fq')):
			tmplfqNumOfLines = getFastQLength(tmplFq)
		else:
			tmplfqNumOfLines = getFastQGzLength(tmplFq)
		
		# Get the number of lines for the vase created fastq file
		vasefqNumOfLines = 0
		if(vaseFq.endswith('fastq') or vaseFq.endswith('.fq')):
			vasefqNumOfLines = getFastQLength(vaseFq)
		else:
			vasefqNumOfLines = getFastQGzLength(vaseFq)
		
		# Calculate the control length the VaSe produced FastQ should have
		controlVaseLength = ((tmplfqNumOfLines - (numOfAcceptorReads*4)) + (numOfDonorReads*4))	# Calculate what the length of the VaSe FastQ should be, multiply with 4 as each read has four lines.
		
		# Check whether the the calculated length (controlVaseLength) and the actual length of the vase fastq are the same. Also check calculated fq length + 1 in case the actual files have an empty line at the end of the file.
		if(controleVaseLength==vasefqNumOfLines or (controleVaseLength+1)==vasefqNumOfLines):
			return True
		return False
	
	
	# Returns the number of lines in a plain text fastq file.
	# Thanks to SilentGhost at https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
	def getFastQLength(fqFileLoc):
		with open(fqFileLoc, 'r') as fqFile:
			for i, l in enumerate(fqFile):
				pass
		return i + 1
	
	
	# Returns the number of lines in a gzipped fastq file.
	def getFastQGzLength(fqFileLoc):
		with gzip.open(fqFileLoc, 'rt') as fqFile:
			for i, l in enumerate(fqFile):
				pass
		return i + 1
