import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class VaSeCompareVariantContext:
    def __init__(self):
        self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
    
    # Run all analysis steps of this util
    def main(self, vcFile1Loc, vcfFile2Loc):
        self.vaseUtilLogger.info("Run VaSe util VaSeCompareVariantContexts")
        varconFile1 = VariantContextFile(vcFile1Loc)
        varconFile2 = VariantContextFile(vcFile2Loc)
        
        # Compare the two varcon files (Maybe add a method to do so in the VariantContextFile itself)
        varconFile1.compare(varconFile2)
        self.vaseUtilLogger.info("Finished running VaSeCompareVariantContexts")
