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
        self.getUnMappedInfo(unmappedData)
        self.vaseutillogger.info("Finished running VaSe util UnmappedMateInfo")

    # Reads the unmapped mate read file
    def read_unmapped_mate_reads(self, unmappedMateLoc, sampleFilter):
        unmappedData = {}
        try:
            with open(unmappedMateLoc, 'r') as urmFile:
                next(urmFile)    # Skip the header line
                
                for fileLine in urmFile:
                    fileLine = fileLine.strip()
                    fileLineData = fileLine.split("\t")
                    unmappedData[fileLineData[0]] = fileLineData[1].split(" ; ")
        except IOError as ioe:
            self.vaseutillogger.critical("Could not read file containing read identifiers with unmapped mates")
        return unmappedData

    # Displays the read info for reads with unmapped mates
    def get_unmapped_info(self):
        for sampleId, readMateInfo
