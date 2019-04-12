import logging

class VaSeUtilHelper:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Returns whether something is in the filter or not
	def passesFilter(self, valToCheck, filterList):
		if(filterList is not None):
			if(valToCheck in filterList):
				return True
			return False
		return True
	
	
	# Reads the file with a list of used donor VCF/BAM files
	def readDonorListFile(self, dListFile, sampleFilter=None):
		donorFiles = {}
		try:
			with open(dListFile, 'r') as dlFile:
				next(dlFile)	# Skip the header line
				for fileLine in dlFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					# Check if the entry is in the set sample filter
					if(self.passesFilter(fileLineData[0], sampleFilter)):
						donorFiles[fileLineData[0]] = fileLineData[1]
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read donor list file")
		return donorFiles
	
	
	# Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck	').
	def readABamReadsList_noFilter(self, acceptorReadFile):
		acceptorReads = []
		with open(acceptorReadFile, 'r') as arFile:
			next(arFile)	# Skip the header line
			for fileLine in arFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split("\t")
				acceptorReads.extend(fileLineData[1:])
		return acceptorReads
	
	
	# Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck	').
	def readDBamReadsList_noFilter(self, donorReadFile):
		donorReads = []
		with open(donorReadFile, 'r') as arFile:
			next(arFile)	# Skip the header line
			for fileLine in arFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split("\t")
				donorReads.extend(fileLineData[2:])
		return donorReads
	
	
	# Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck	').
	def readDBamReadsListByVarcon_noFilter(self, acceptorReadFile):
		donorReads = []
		with open(acceptorReadFile, 'r') as arFile:
			next(arFile)	# Skip the header line
			for fileLine in arFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split("\t")
				donorReads.extend(fileLineData[2:])
		return donorReads
