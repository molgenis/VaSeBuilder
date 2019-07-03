import logging
from VaSeUtilHelper import VaSeUtilHelper

class LogInfo:
    def __init__(self, vaseuhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseuhelper

    # Runs the loginfo program and prints info from the log file.
    def main(self, vaselogloc, logfilter):
        self.vaseutillogger.info("Running VaSe util LogInfo")
        self.process_logfile(vaselogloc, logfilter)
        self.vaseutillogger.info("Finished running VaSe util LogInfo")

    # Processes the log file and prints all lines or lines satisfying the filter
    def process_logfile(self, vaselogloc, logfilter=None):
        try:
            with open(vaselogloc, 'r') as vaselogfile:
                for fileline in vaselogfile:
                    fileline_elements = fileline.split("\t")
                    if self.vuh.passes_filter(fileline_elements[2], logfilter):
                        print(fileline.strip())
        except IOError as ioe:
            self.vaseutillogger.warning("Could not open log file " + str(ioe.filename))
