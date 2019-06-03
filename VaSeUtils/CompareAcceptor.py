import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareAcceptor:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
    
    # Compares the acceptor reads for two VaSe outputs
    def main(self, varconfileloc1, varconfileloc2):
        self.vaseutillogger.info("Running VaSE util CompareAcceptor")
        varconFile1 = VariantContextFile(varconfileloc1)
        varconFile2 = VariantContextFile(varconfileloc2)
        self.performComparison(varconFile1, varconFile2)
        self.vaseutillogger.info("Finished running VaSe util CompareAcceptor")

    # Performs the comparison per variant context
    def performComparison(self, varconFile1, varconFile2):
        varconAcceptorReads1 = varconFile1.getAllAcceptorReadIdsByVarcon()
        varconAcceptorReads2 = varconFile2.getAllAcceptorReadIdsByVarcon()
        
        varContextsNotInOther = []
        if(len(varconAcceptorReads2) > len(varconAcceptorReads1)):
            self.vaseutillogger.info("Comparing acceptor reads from the second varcon file to the first")
            self.compareAcceptorReads()
        else:
            self.vaseutillogger.info("Comparing acceptor reads from the first varcon file to the second")

    # Performs the comparison per variant context
    def perform_comparison(self, varconfile1, varconfile2):
        acceptorreads1 = varconfile1.get_all_variant_context_acceptor_reads()
        acceptorreads2 = varconfile2.get_all_variant_context_acceptor_reads()

        print(f"The first variant context file contains {len(acceptorreads1)} acceptor reads")
        print(f"The second variant context file contains {len(acceptorreads2)} acceptor reads")

        if len(acceptorreads1) > 0 and len(acceptorreads2) > 0:
            self.vaseutillogger.info("Comparing variant context acceptor reads between two variant context files")
            self.compare_acceptor_reads()

    # Iterates over the largest map with reads per variant context and compares to the smaller one
    def compareAcceptorReads(self, largerReadData, smallerReadData):
        contextsNotInOther = []
        for varconId, acceptorReads in largerReadData.items():
            if():
                
            else:
                contextsNotInOther.append

    # Compares the acceptor reads
    def compare_acceptor_reads(self, areadlist1, areadlist2):
