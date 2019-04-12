import logging
import pysam

from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class DonorReadInfo:
	def __init__(self, vaseuhelper):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = vaseuhelper
	
	
	# Performs the main analysis
	def main(self, donorBamListFile, donorbreadFile, vcFile, sampleFilter=None, varconFilter=None):
		self.vaseUtilLogger.inf("Running VaSe util DonorReadInfo")
		dbamFiles = self.readDonorBamListFile(donorBamListFile, sampleFilter)
		dbreads = self.readDonorBreadFile(donorbreadFile, sampleFilter, varconFilter)
		varconFile = VariantContextFile(vcFile, sampleFilter, varconFilter)
		self.getDonorReadInfo(dbamFiles, dbreads, varconFile)
		self.vaseUtilLogger.info("Finished running VaSe util DonorReadInfo")
	
	
	# Reads the list of used donor BAM files
	def readDonorBamListFile(self, dbamListFile, sampleFilter):
		donorBams = {}
		try:
			with open(dbamListFile, 'r') as dblFile:
				next(dblFile)	# Skip the header line
				for fileLine in dblFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					# Check if the entry is in the set sample filter
					if(self.passesFilter(fileLineData[0], sampleFilter)):
						donorBams[fileLineData[0]] = fileLineData[1]
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read")
		return donorBams
	
	
	# Reads the donorbread file
	def readDonorBreadFile(self, donorbreadFile, sampleFilter, varconFilter):
		donorBreads = {}
		try:
			with open(donorbreadFile, 'r') as dbrFile:
				next(dbrFile)	# Skip the header line
				for fileLine in dbrFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.vuh.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.vuh.passesFilter(fileLineData[0], varconFilter)
					
					# Add the read data to the map
					if(samplePass and varconPass):
						if(fileLineData[1] not in donorBreads):
							donorBreads[1] = {}
						donorBreads[fileLineData[1]][fileLineData[0]] = fileLineData[2].split(" ; ")
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read the donorbread file.")
		return donorBreads
	
	
	# Obtains the read info for the donor BAM reads satisfying the set sample and variant context filters
	def getDonorReadInfo(self, donorBreads, donorBams, varconFile):
		for sampleId, varconReads in donorBreads.items():
			self.vaseUtilLogger.info("SAMPLE: " +str(sampleId))
			for varconId, dbreads in varconReads.items():
				searchChrom = varconFile.getVariantContextChrom(varconId)
				searchStart = varconFile.getVariantContextStart(varconId)
				searchStop = varconFile.getVariantContextEnd(varconId)
				
				# If all three required BAM searching parameters are valid start searching for the reads
				if(searchChrom and searchStart and searchEnd):
					try:
						dBamFile = pysam.AlignmentFile(donorBams[sampleId])
						self.vaseUtilLogger.info("Read info for variant context: " +str(varconId))
						for bread in dBamFile.fetch(searchChrom, searchStart, searchStop):
							if(bread.query_name in dbreads):
								self.vaseUtilLogger.info(bread.to_string())
						dBamFile.close()
					except IOError as ioe:
						self.vaseUtilLogger.warning("Could not read BAM file")