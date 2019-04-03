import logging
import pysam

from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class DonorReadInfo:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Performs the main analysis
	def main(self, donorBamListFile, donorbreadFile, vcFile, sampleFilter=None, varconFilter=None):
		dbamFiles = self.readDonorBamListFile(donorBamListFile, sampleFilter)
		dbreads = self.readDonorBreadFile(donorbreadFile, sampleFilter, varconFilter)
		varconFile = VariantContextFile(vcFile, sampleFilter, varconFilter)
		self.getDonorReadInfo(dbamFiles, dbreads, varconFile)
	
	
	# Reads the list 
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
			return donorBams
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read")
	
	
	# Reads the donorbread file
	def readDonorBreadFile(self, donorbreadFile, sampleFilter, varconFilter):
		donorBreads = {}
		try:
			with open(donorbreadFile, 'r') as dbrFile:
				next(dbrFile)	# Skip the header line
				for fileLine in dbrFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					
					samplePass = self.passesFilter(fileLineData[1], sampleFilter)
					varconPass = self.passesFilter(fileLineData[0], varconFilter)
					
					# Add the read data to the map
					if(samplePass and varconPass):
						if(fileLineData[1] not in donorBreads):
							donorBreads[1] = {}
						donorBreads[fileLineData[1]][fileLineData[0]] = fileLineData[2].split(" ; ")
			return donorBreads
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read the donorbread file.")
	
	
	# Returns whether something is in the filter
	def passesFilter(self, valToCheck, filterList):
		if(filterList is not None):
			if(valToCheck in filterList):
				return True
			return False
		return True
	
	
	# Obtains the read info for the donor BAM reads satisffying the set sample and variant context filters
	def getDonorReadInfo(self, donorBams, donorBreads):
		for sampleId, dBams in donorBams.items():
			for dBam in dBams:
				try:
					dbamFile = pysam.AlignmentFile(dBam)
					
				except IOError as ioe:
					self.vaseUtilLogger.warning("")
	
	
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
								self.vaseUtilLogger.info(bread.query_name+"\t"+bread)
								#Print info for each bamread per line
						dBamFile.close()
					except IOError as ioe:
						self.vaseUtilLogger.warning("Could not read BAM file")
