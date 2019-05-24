import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class UnmappedMateInfo:
    def __init__(self, vaseuhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Performs the unmapped mate read info analysis
    def main(self, unmappedmateloc, acceptorBamLoc, donorfilesloc, samplefiltr=None, varconFilter=None):
        self.vaseutillogger.info("Running VaSe util UnmapppedMateInfo")
        varconFile = VariantContextFile()
        
        unmappedData = self.read_unmapped_mate_reads(unmappedmateloc)
        donorBamFiles = self.dri.readDonorBamListFile(donorfilesloc, samplefiltr)
        self.get_unmapped_info(unmappedData)
        self.vaseutillogger.info("Finished running VaSe util UnmappedMateInfo")

    # Reads the unmapped mate read file
    def read_unmapped_mate_reads(self, unmappedmateloc, sampleFilter):
        unmappeddata = {}
        try:
            with open(unmappedmateloc, 'r') as urmfile:
                next(urmfile)    # Skip the header line
                
                for fileline in urmfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    unmappeddata[filelinedata[0]] = filelinedata[1].split(" ; ")
        except IOError as ioe:
            self.vaseutillogger.critical("Could not read file containing read identifiers with unmapped mates")
        return unmappeddata

    # Displays the read info for reads with unmapped mates
    def get_unmapped_info(self):
        for sampleId, readMateInfo
