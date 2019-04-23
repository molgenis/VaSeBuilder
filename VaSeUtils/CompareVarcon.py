import logging
from VariantContext impoort VariantContext
from VariantContextFile import VariantContextFile

class CompareVarcon
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Performs the main analysis
	def main(self, varconFileLoc1, varconFileLoc2):
		seflf.vaseUtilLogger.info("Running VaSe util CompareVarcon")
		varconFile1 = VariantContextFile(varconFileLoc1)
		varconFile2 = VariantContextFile(varconFileLoc2)
		self.compareVarconFiles(varconFile1, varconFile2)
		self.vaseUtilLogger.info("Finished running VaSeUtil CompareVarcon")
	
	
	# Performs the comparison between two variant context files.
	def compareVarconFiles(self, varconFile1, varconFile2):
		self.displayGeneralInfo(varconFile1, varconFile2)
	
	
	# Displays general info about the two variant contexts files
	def displayGeneralInfo(self, varconFile1, varconFile2):
		print("[-Variant Context General Info-]")
		print("\tVarcon1\tVarcon2")
		print("#Contexts:\t" +str(varconFile1.getNumberOfVariantContexts())+ "\t" +str(varconFile2.getNumberOfVariantContexts()))
		print("Union:")
		print("Intersect: ")
		print("")
