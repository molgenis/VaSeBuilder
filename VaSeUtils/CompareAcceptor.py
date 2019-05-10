import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareAcceptor:
	def __init__(self):
		self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
	
	# Compares the acceptor reads for two VaSe outputs
	def main(self, varconFileLoc1, varconFileLoc2):
		self.vaseUtilLogger.info("Running VaSE util CompareAcceptor")
		varconFile1 = VariantContextFile(varconFileLoc1)
		varconFile2 = VariantContextFile(varconFileLoc2)
		self.performComparison(varconFile1, varconFile2)
		self.vaseUtilLogger.info("Finished running VaSe util CompareAcceptor")
	
	
	# Performs the comparison per variant context
	def performComparison(self, varconFile1, varconFile2):
		varconAcceptorReads1 = varconFile1.getAllAcceptorReadIdsByVarcon()
		varconAcceptorReads2 = varconFile2.getAllAcceptorReadIdsByVarcon()
		
		varContextsNotInOther = []
		if(len(varconAcceptorReads2) > len(varconAcceptorReads1)):
			self.vaseUtilLogger.info("Comparing acceptor reads from the second varcon file to the first")
			self.compareAcceptorReads()
		else:
			self.vaseUtilLogger.info("Comparing acceptor reads from the first varcon file to the second")
	
	
	# Iterates over the largest map with reads per variant context and compares to the smaller one
	def compareAcceptorReads(self, largerReadData, smallerReadData):
		contextsNotInOther = []
		for varconId, acceptorReads in largerReadData.items():
			if():
				
			else:
				contextsNotInOther.append
