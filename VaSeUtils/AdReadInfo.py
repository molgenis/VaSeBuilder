import logging
import pysam
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class AdReadInfo:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaseUtil_Logger")
	
	
	# Performs the main Acceptor/Donor Read info work.
	def main(self, varconFileLoc, adReadFile, bamFile, sampleFilter=None, varconFilter=None, chromFilter=None, posFilter=None):
		vcFile = VariantContextFile(varconFileLoc, sampleFilter, varconFilter, chromFilter, posFilter)
		bamReadList = self.readAdReadsList(adReadFile)
		self.getReadInfo(vcFile.getVariantContextsBySample(), bamReadList, acceptm)
	
	
	# Reads the file containing acceptor/donor BAM reads.
	def readAdReadsList(self, adReadFile):
		adReads = []
		with open(adReadFile, 'r') as arFile:
			next(arFile)	# Skip the header line
			for fileLine in arFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split("\t")
				adReads.extend(fileLineData[1:])
		return adReads
	
	
	# Reads the file containing acceptor
	def readAdReadsList(self, adReadFile, sampleFilter, varconFilter, chromFilter, posFilter):
	
	
	# Reads the file containing the used donor BAM files.
	def readDonorBamListFile(self, dbamListFile):
		donorBams = {}
		try:
			with open(dbamListFile, 'r') as dbReadFile:
				next(dbReadFile)
				for fileLine in dbReadFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
						
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read donorbread file.")
	
	
	# Obtains the BAM information about the read
	def getReadInfo(self, varconData, readsList, bamFile):
		
