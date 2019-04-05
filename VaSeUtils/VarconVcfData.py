import logging
import pysam
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class VarconVcfData:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vaseUtilLogger.info("Running VaSe util VarconData")
	
	
	# Runs all the analysis steps
	def main(self, dvcfListFileLoc, varconFileLoc, sampeFilter=None, varconFilter=None, chromFilter=None):
		self.readUsedDonorListFile(dvcfListFile, sampleFilter)
		varconFile = VariantContextFile(varconFileLoc, sampleFilter, varconFilter, chromFilter)
		self.getVarconVcfData()
	
	
	# Reads the donor list file and saves it into a 
	def readUsedDonorListFile(dListFileLoc, sampleFilter=None):
		donorList = {}
		try:
			with open(dListFileLoc, 'r') as dlistFile:
				next(dlistFile)	# Skip the header line
				for fileLine in dlistFile:
					fileLine = fileLine.strip()
					fileLineData = fileLine.split()
					
					if(self.passesFilter(fileLineData[0])):
						donorList[fileLineData[0]] = fileLineData[1:]
		except IOError as ioe:
			print("Could not read donor list file.")
			exit()
		return donorList
	
	
	# Displays the VCF data for selected variant contexts (all if no filters were set)
	def getVarconVcfData(dListFileLoc, varconFileLoc, sampleFilter=None, chromFilter=None, posFilter=None):
		donorList = readUsedDonorListFile(dListFileLoc)
		varconFile = VariantContextFile(varconFileLoc, sampleFilter, chromFilter, posFilter)
		variantContexts = varconFile.getVariantContexts()
		
		for dvcfFileSample in donorList:
			processVcfFile(dvcfFileSample, donorList[dvcfFileSample])
	
	
	# Processes a single VCF file.
	def processVcfFile(sampleId, varconData, vcfFileLoc):
		try:
			vcfFile = pysam.VariantFile(vcfFileLoc)
			for varcon in varconData:
				vcfFile.fetch(varcon.getVariantContextChrom(), )
			vcfFile.close()
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read VCF file for sample " +str(vcfFileLoc))
	
	
	# Returns whether something is in the filter
	def passesFilter(self, valToCheck, filterList):
		if(filterList is not None):
			if(valToCheck in filterList):
				return True
			return False
		return True
