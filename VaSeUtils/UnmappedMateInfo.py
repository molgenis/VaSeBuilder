import logging
from AcceptorReadInfo import AccceptorReadInfo
from DonorReadInfo import DonorReadInfo

class UnmappedMateInfo:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.ari = AcceptorReadInfo()
		sel.dri = DonorReadInfo()
	
	
	# Performs the unmapped mate read info analysis
	def main(self, unmappedMateLoc, acceptorBamLoc, donorFilesLoc, sampleFilter=None):
		self.vaseUtilLogger.info("Running VaSe util UnmapppedMateInfo")
		unmappedData = self.readUnmappedMateReads(unmappedMateLoc)
		donorBamFiles = self.dri.readDonorBamListFile(donorFilesLoc, sampleFilter)
		self.getUnMappedInfo(unmappedData)
		self.vaseUtilLogger.info("Finished running VaSe util UnmappedMateInfo")
	
	
	# Reads the unmapped mate read file
	def readUnmappedMateReads(self, unmappedMateLoc, sampleFilter):
		unmappedData = {}
		try:
			with open(unmappedMateLoc, 'r') as urmFile:
				next(urmFile)	# Skip the header line
				
				for fileLine in urmFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					unmappedData[fileLineData[0]] = fileLineData[1].split(" ; ")
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read file containing read identifiers with unmapped mates")
		return unmappedData
	
	
	# Displays the read info for reads with unmapped mates
	def getUnmappedInfo(self):
		for sampleId, readMateInfo
