import logging
import pysam

from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class AcceptorReadInfo:
	def __init__(self, vaseuhelper):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = vaseuhelper
	
	
	# Performs all the analysis steps
	def main(self, acceptorBamFile, acceptorbreadFile, vcFileLoc, sampleFilter=None, varconFilter=None):
		self.vaseUtilLogger.info("Running VaSe util AcceptorReadInfo")
		abreads = self.readAcceptorBreadFile(acceptorbreadFile, sampleFilter, varconFilter)
		varconFile = VariantContextFile(vcFileLoc, sampleFilter, varconFilter)
		self.getDonorReadInfo(abreads, acceptorBamFile, varconFile)
		self.vaseUtilLogger.info("Finished running VaSe util AcceptorReadInfo")
	
	
	# Reads the acceptorbread file 
	def readAcceptorBreadFile(self, acceptorbreadFile, sampleFilter, varconFilter):
		acceptorBreads = {}
		try:
			with open(acceptorbreadFile, 'r') as abrFile:
				next(abrFile)	# Skip the header line
				for fileLine in dbrFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.vuh.passesFilter(fileLineData[1], None)
					varconPass = self.vuh.passesFilter(fileLineData[0], varconFilter)
					
					# Add the read data to the map
					if(samplePass and varconPass):
						acceptorBreads[fileLineData[0]] = fileLineData[1].split(" ; ")
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read the acceptorbread file.")
		return acceptorBreads
	
	
	# Obtains the read info for the selected acceptor reads (all if no filters were set all reads will be used)
	def getAcceptorReadInfo(self, acceptorBreads, acceptorBam, varconFile):
		try:
			aBamFile = pysam.AlignmentFile(acceptorBam)
			for varconId, varconReads in acceptorBreads.items():
				searchChrom = varconFile.getVariantContextChrom(varconId)
				searchStart = varconFile.getVariantContextStart(varconId)
				searchStop = varconFile.getVariantContextEnd(varconId)
				
				if(searchChrom and searchStart and searchStop):
					for abread in aBamFile.fetch(searchChrom, searchStart, searchStop):
						if(abread.query_name in varconReads):
							self.vaseUtilLogger.info(bread.to_string())
			aBamFile.close()
		except IOError as ioe:
