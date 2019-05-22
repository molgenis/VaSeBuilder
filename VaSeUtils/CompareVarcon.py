import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareVarcon:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Performs the main analysis
	def main(self, varconFileLoc1, varconFileLoc2):
		self.vaseUtilLogger.info("Running VaSe util CompareVarcon")
		varconFile1 = VariantContextFile(varconFileLoc1)
		varconFile2 = VariantContextFile(varconFileLoc2)
		self.compareVarconFiles(varconFile1, varconFile2)
		self.vaseUtilLogger.info("Finished running VaSeUtil CompareVarcon")
	
	
	# Performs the comparison between two variant context files.
	def compareVarconFiles(self, varconFile1, varconFile2):
		self.displayGeneralInfo(varconFile1, varconFile2)
		self.displaySharedVarconDifferences(varconFile1, varconFile2)
	
	
	# Displays general info about the two variant contexts files
	def displayGeneralInfo(self, varconFile1, varconFile2):
		print("==================================================")
		print("[-Variant Context General Info-]")
		print("\tVarcon1\tVarcon2")
		print("#Contexts:\t" +str(varconFile1.getNumberOfVariantContexts())+ "\t" +str(varconFile2.getNumberOfVariantContexts()))
		print("#Difference:\t" + str(varconFile1.get_variant_contexts_difference(varconFile2)) + "\t" + str(varconFile2.get_variant_contexts_difference(varconFile1)))
		print("--------------------")
		print("#Union: " + str(varconFile1.get_variant_contexts_union(varconFile2)))
		print("#Intersect: " + str(varconFile1.get_variant_contexts_intersect(varconFile2)))
		print("#Symmetric Difference: " + str(varconFile1.get_variant_contexts_symmetric_difference(varconFile2)))
	
	
	# Displays which variant contexts are in both but differ (context start/end, etc)
	def displaySharedVarconDifferences(self, varconFile1, varconFile2):
		variantContexts1 = varconFile1.getVariantContextsById()
		variantContexts2 = varconFile2.getVariantContextsById()
		sharedVarconList = varconFile1.get_variant_contexts_intersect(varconFile2)
		sharedVarconDiffs = varconFile1.compare(varconFile2, sharedVarconList)
		diffNumberMap = varconFile1.getVariantContextFields()
		
		# Iterate over the shared variant contexts differences
		print("==================================================")
		print("[-Differing variant contexts-]")
		for varconId, varconDiffs in sharedVarconDiffs.items():
			msg= str(varconId)+": "
			for diffield, diffvalues in varconDiffs:
				msg += str(diffNumberMap[diffield])+ " (" +str(diffvalues[0])+ "|" +str(diffvalues[1])+ "), "
			print(msg)
