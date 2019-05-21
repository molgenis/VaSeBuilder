import logging
import pysam
from VaSeUtilHelper import VaSeUtilHelper
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class VarconVcfData:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
		self.vuh = VaSeUtilHelper()
	
	
	# Runs all the analysis steps
	def main(self, dvcfListFileLoc, varconFileLoc, sampeFilter=None, varconFilter=None, chromFilter=None):
		self.vaseUtilLogger.info("Running VaSe util VarconData")
		self.getVarconVcfData(dvcfListFileLoc, varconFileLoc, sampeFilter, varconFilter, chromFilter)
		self.vaseUtilLogger.info("Finished running VaSe util VarconData")
	
	
	# Displays the VCF data for selected variant contexts (all if no filters were set)
	def getVarconVcfData(dListFileLoc, varconFileLoc, sampleFilter=None, varconFilter=None, chromFilter=None):
		donorList = self.vuh.readUsedDonorListFile(dListFileLoc, sampleFilter)
		varconFile = VariantContextFile(varconFileLoc, sampleFilter, varconFilter, chromFilter)
		variantContexts = varconFile.getVariantContexts()
		
		for dvcfFileSample in donorList:
			print("Varcon\tChrom\tPos\tRef")
			processVcfFile(dvcfFileSample, variantContexts, donorList[dvcfFileSample])
	
	
	# Processes a single VCF file.
	def processVcfFile(sampleId, varconData, vcfFileLoc):
		try:
			vcfFile = pysam.VariantFile(vcfFileLoc)
			for varcon in varconData:
				for dvcfVar in vcfFile.fetch(varcon.get_variant_context_chrom(), varcon.get_variant_context_origin(), varcon.get_variant_context_origin()):
					if(varcon.get_variant_context_origin() == dvcfVar.pos):
						print(str(varcon.get_variant_context_id()) + "\t" + str(dvcfVar.chrom) + "\t" + str(dvcfVar.pos) + "\t" + str(dvcfVar.ref))
			vcfFile.close()
		except IOError as ioe:
			self.vaseUtilLogger.critical("Could not read VCF file for sample " +str(vcfFileLoc))
