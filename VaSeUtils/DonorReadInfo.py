import logging
import pysam

# Import required VaSeUtil classes
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class DonorReadInfo:
	def __init__(self, vaseuhelper):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = vaseuhelper
	
	# Performs the main analysis
	def main(self, donorBamListFile, donorbreadFile, vcFile, sampleFilter=None, varconFilter=None, readIdFilter=None):
		self.vaseUtilLogger.info("Running VaSe util DonorReadInfo")
		dbamFiles = self.vuh.readDonorListFile(donorBamListFile, sampleFilter)
		varconFile = VariantContextFile(vcFile, sampleFilter, varconFilter)
		dbreads = varconFile.getAllDonorReadIdsByVarcon()
		self.getDonorReadInfo(dbamFiles, dbreads, varconFile, readIfFilter)
		self.vaseUtilLogger.info("Finished running VaSe util DonorReadInfo")
	
	
	# Obtains the read info for the donor BAM reads satisfying the set sample and variant context filters
	def getDonorReadInfo(self, donorBreads, donorBams, varconFile, readIdFilter):
		for sampleId, varconReads in donorBreads.items():
			print("SAMPLE: " +str(sampleId))
			for varconId, dbreads in varconReads.items():
				searchChrom = varconFile.getVariantContextChrom(varconId)
				searchStart = varconFile.getVariantContextStart(varconId)
				searchStop = varconFile.getVariantContextEnd(varconId)
				
				# If all three required BAM searching parameters are valid start searching for the reads
				if(searchChrom and searchStart and searchEnd):
					try:
						dBamFile = pysam.AlignmentFile(donorBams[sampleId])
						print("Read info for variant context: " +str(varconId))
						for bread in dBamFile.fetch(searchChrom, searchStart, searchStop):
							if(bread.query_name in dbreads):
								self.vaseUtilLogger.info(bread.to_string())
						dBamFile.close()
					except IOError as ioe:
						self.vaseUtilLogger.warning("Could not read BAM file")
