import logging

class CompareVarconStats:
    def __init__(self, varconStatsLoc):
        self.varconStatsData = self.readVarconStatsFile(varconStatsLoc)
    
    
    # Read the Variant Context status file
    def readVarconStatsFile(self, fileLoc):
        statsData = {}
        try:
            with open(fileLoc, 'r') as varconStatsFile:
                next(varconStatsFile)    # Skip the header line
                
                for fileLine in varconStatsFile:
                    fileLine = fileLine.strip()
                    fileLineData = fileLine.split("\t")
        except IOError as ioe:
            self.vaseUtilLogger.warning("Could not read " +str())
        return statsData
