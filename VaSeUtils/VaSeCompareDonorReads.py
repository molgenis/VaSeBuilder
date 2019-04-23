import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class VaSeCompareDonorReads:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	
	# Compares the list of donor reads from two donorbread.txt files.
	def main(self, varconLoc1, varconLoc2):
		self.vaseUtilLogger.info("Running VaSe util VaSeCompareDonorReads")
		varconFile1 = VariantContextFile(varconLoc1)
		varconFile2 = VariantContextFile(varconLoc2)
		dReadsList1 = varconFile.getAllDonorReadIdsByVarcon()
		dReadsList2 = varconFile.getAllDonorReadIdsByVarcon()
		
		# Check which read list is the largest.
		if(len(dReadList1) > len(dReadList2)):
			self.compareDonorReads(dReadList1, dReadList2)
		else:
			self.compareDonorReads(dReadList2, dReadList1)
		self.vaseUtilLogger.info("Finished running util VaSeCompareDonorReads")
	
	
	# Compares tow lists of donor reads from donorbread files from two runs
	def compareDonorReads(self, dreadsL, dreadsS):
		for varconId, dReads in dreadsL.items():
			if(varconId in dreads2):
				vcReadDiffs = self.compareReadLists(dreadsL[varconId], dreadsS[varconId])
				if(len(vcReadDiffs)>0):
					self.vaseUtilLogger.info("Variant context " +str(varconId)+ " differs by " +str()+ " reads.")
					self.vaseUtilLogger.info("Differing reads: " +" ; ".join(vcReadDiffs))
			else:
				self.vaseUtilLogger.info("No donor reads for variant context " +str(varconId)+ " in the second donor read file")
	
	
	# Thanks to Roman Bodnarchuk at https://stackoverflow.com/questions/6486450/python-compute-list-difference
	def getReadListDifferences(readList1, readList2):
		readList2 = set(readList2)
		return [readId for readId in readList1 if readId not in readList2]
