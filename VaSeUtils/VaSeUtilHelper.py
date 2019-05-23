import logging

class VaSeUtilHelper:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Returns whether something is in the filter or not
    def passes_filter(self, valtocheck, filterlist):
        if filterlist is not None:
            if valtocheck in filterlist:
                return True
            return False
        return True

    # Reads the file with a list of used donor VCF/BAM files
    def read_donor_list_file(self, dlistfile, samplefilter=None):
        donorfiles = {}
        try:
            with open(dlistfile, 'r') as dlfile:
                next(dlfile)    # Skip the header line
                for fileline in dlfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    
                    # Check if the entry is in the set sample filter
                    if self.passes_filter(filelinedata[0], samplefilter):
                        donorfiles[filelinedata[0]] = filelinedata[1]
        except IOError as ioe:
            self.vaseutillogger.critical("Could not read donor list file")
        return donorfiles

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck    ').
    def read_abamreads_list_nofilter(self, acceptorReadFile):
        acceptorReads = []
        with open(acceptorReadFile, 'r') as arFile:
            next(arFile)    # Skip the header line
            for fileLine in arFile:
                fileLine = fileLine.strip()
                fileLineData = fileLine.split("\t")
                acceptorReads.extend(fileLineData[1:])
        return acceptorReads

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck    ').
    def readDBamReadsList_noFilter(self, donorReadFile):
        donorReads = []
        with open(donorReadFile, 'r') as arFile:
            next(arFile)    # Skip the header line
            for fileLine in arFile:
                fileLine = fileLine.strip()
                fileLineData = fileLine.split("\t")
                donorReads.extend(fileLineData[2:])
        return donorReads

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck    ').
    def readDBamReadsListByVarcon_noFilter(self, acceptorReadFile):
        donorReads = []
        with open(acceptorReadFile, 'r') as arFile:
            next(arFile)    # Skip the header line
            for fileLine in arFile:
                fileLine = fileLine.strip()
                fileLineData = fileLine.split("\t")
                donorReads.extend(fileLineData[2:])
        return donorReads
