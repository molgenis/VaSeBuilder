import os
import logging

class CheckDonorFilesExist:
    def __init__(self, vaseuhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseuhelper

    def main(self, donorlistfile, samplefilter=None):
        self.vaseutillogger.info("Running VaS util CheckDonorFilesExist")
        try:
            with open(donorlistfile, "r") as dlfile:
                next(dlfile)    # Skip the header line of the file
                for fileline in dlfile:
                    filelinedata = fileline.strip().split("\t")
                    donorfiles = filelinedata[1].split(",")

                    # Check whether the used donor files exist
                    if self.vuh.passes_filter(filelinedata[0], samplefilter):
                        for dfile in donorfiles:
                            if not os.path.isfile(dfile):
                                print(f"Donor file {dfile} could not be found. Maybe it has been moved, renamed or"
                                      "deleted?")
        except IOError:
            self.vaseutillogger.warning(f"Could not open {donorlistfile}")
        self.vaseutillogger.info("Finished running VaSe util CheckDonorFilesExist")
