import logging
import pysam

# Import required VaSeUtils classes
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class AcceptorReadInfo:
    def __init__(self, vaseuhelper):
        self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseuhelper
    
    
    # Performs all the analysis steps
    def main(self, acceptorBamFile, vcFileLoc, sampleFilter=None, varconFilter=None, readIdFilter=None):
        self.vaseUtilLogger.info("Running VaSe util AcceptorReadInfo")
        varconFile = VariantContextFile(vcFileLoc, sampleFilter, varconFilter)
        abreads = varconFile.getAllAcceptorReadIdsByVarcon()
        self.getAcceptorReadInfo(abreads, acceptorBamFile, varconFile)
        self.vaseUtilLogger.info("Finished running VaSe util AcceptorReadInfo")
    
    
    # Obtains the read info for the selected acceptor reads (all if no filters were set all reads will be used)
    def getAcceptorReadInfo(self, acceptorBreads, acceptorBam, varconFile):
        try:
            aBamFile = pysam.AlignmentFile(acceptorBam)
            for varconId, varconReads in acceptorBreads.items():
                searchChrom = varconFile.get_variant_context_chrom(varconId)
                searchStart = varconFile.get_variant_context_start(varconId)
                searchStop = varconFile.get_variant_context_end(varconId)
                
                if(searchChrom and searchStart and searchStop):
                    print("Variant Context: " + str(varconId) + " ;; from sample: " + str(varconFile.get_sample_id(varconId)))
                    for abread in aBamFile.fetch(searchChrom, searchStart, searchStop):
                        if(abread.query_name in varconReads):
                            if(self.vuh.passes_filter(abread.query_name, readIdFilter)):
                                print(bread.to_string())
                    print("")
            aBamFile.close()
        except IOError as ioe:
