import logging
import pysam

class AcceptorReadInfo:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Performs all the analysis steps
	def main(self, acceptorBamFile, acceptorbreadFile, vcFileLoc, varconFilter=None):
		self.vaseUtilLogger.info("Running VaSe util AcceptorReadInfo")
		abreads = self.readDonorBreadFile(donorbreadFile, sampleFilter, varconFilter)
		varconFile = VariantContextFile(vcFileLoc, sampleFilter, varconFilter)
		self.getDonorReadInfo(acceptorBamFile, abreads, varconFile)
		self.vaseUtilLogger.info("Finished running VaSe util AcceptorReadInfo")
	
	
	# Reads the acceptorbread file 
	def readAcceptorBreadFile(self, donorbreadFile, sampleFilter, varconFilter):
		acceptorBreads = {}
		try:
			with open(donorbreadFile, 'r') as dbrFile:
				next(dbrFile)	# Skip the header line
				for fileLine in dbrFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.passesFilter(fileLineData[1], None)
					varconPass = self.passesFilter(fileLineData[0], varconFilter)
					
					# Add the read data to the map
					if(samplePass and varconPass):
						acceptorBreads[fileLineData[0]] = fileLineData[1].split(" ; ")
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read the acceptorbread file.")
		return acceptorBreads
	
	
	# Returns whether something is in the filter or not
	def passesFilter(self, valToCheck, filterList):
		if(filterList is not None):
			if(valToCheck in filterList):
				return True
			return False
		return True
	
	
	# Obtains the read info for the selected acceptor reads (all if no filters were set all reads will be used)
	def getAcceptorReadInfo(self, acceptorBreads, acceptorBam, varconFile):
		aBamFile = pysam.AlignmentFile(acceptorBam)
		for varconId, varconReads in acceptorBreads.items():
			searchChrom = varconFile.getVariantContextChrom(varconId)
			searchStart = varconFile.getVariantContextStart(varconId)
			searchStop = varconFile.getVariantContextEnd(varconId)
			
			if(searchChrom and searchStart and searchStop):
				for abread in aBamFile.fetch(searchChrom, searchStart, searchStop):
					if(abread.query_name in ):
		aBamFile.close()
