import os
import logging

"""Checks whether the used donor files still exist at the recorded location"""


class CheckDonorFilesExist:
    def __init__(self, vaseuhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseuhelper

    # Performs the actual checking if the donor files still exist
    def main(self, donorlistfile, samplefilter=None):
        self.vaseutillogger.info("Running VaS util CheckDonorFilesExist")
        missing_count = 0
        total_count = 0
        try:
            with open(donorlistfile, "r") as dlfile:
                next(dlfile)    # Skip the header line of the file
                for fileline in dlfile:
                    filelinedata = fileline.strip().split("\t")
                    donorfiles = filelinedata[1].split(",")

                    # Check whether the used donor files exist
                    if self.vuh.passes_filter(filelinedata[0], samplefilter):
                        for dfile in donorfiles:
                            total_count += 1
                            if not os.path.isfile(dfile):
                                missing_count += 1
                                print(f"Donor file {dfile} for sample {filelinedata[0]} could not be found. Maybe it "
                                      "has been moved, renamed or deleted?")
            print(f"Found {total_count - missing_count}/{total_count} donor files.")
        except IOError:
            self.vaseutillogger.critical(f"Could not open {donorlistfile}")
        finally:
            self.vaseutillogger.info("Finished running VaSe util CheckDonorFilesExist")
