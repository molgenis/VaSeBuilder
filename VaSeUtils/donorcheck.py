import logging
import gzip

class DonorCheck:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaseUtil_Logger")

    # Performs the main analysis
    def main(self, readslist, vasefq1, vasefq2):
        self.vaseutillogger.info("Running VaSe util DonorCheck")

        # Check if the donor reads have indeed been added to the VaSe R1 FastQ
        self.vaseutillogger.info("Checking the R1 VaSe FastQ file.")
        r1_added = self.get_number_of_reads_added(vasefq1, readslist)
        self.vaseutillogger.info("Added " + str(r1_added) + " of " + str(len(readslist))
                                 + " to the R1 VaSe FastQ file.")

        # Check if the donor reads have indeed been added to the VaSe R2 FastQ
        self.vaseutillogger.info("Checking the R2 VaSe FastQ file.")
        r2_added = self.get_number_of_reads_added(vasefq2, readslist)
        self.vaseutillogger.info("Added " + str(r2_added) + " of " + str(len(readslist))
                                 + " to the R2 VaSe FastQ file.")
        
        if r1_added == len(readslist) and r2_added == len(readslist):
            self.vaseutillogger.info("All donor reads have been added to the VaSe FastQ files")
        else:
            self.vaseutillogger.info("Not all donor reads have been added to the VaSe FastQ files")
        self.vaseutillogger.info("Finished running VaSe util DonorCheck")

    # Returns the number of added reads to R1/R2 VaSe fastq data.
    def get_number_of_reads_added(self, vasefqfiles, readidlist):
        addedreads = 0
        for vfqFile in vasefqfiles:
            addedreads = self.check_donor_reads_added(vfqFile, readidlist, addedreads)
        return addedreads

    # Checks whether the list of donor reads are indeed added to a specified fastq file based on read identifier.
    def check_donor_reads_added(self, gzresultsfile, donorreadlist, addedcount):
        with gzip.open(gzresultsfile, 'rt') as gzfile:
            for fileline in gzfile:
                fileline = fileline.strip()
                if fileline.startswith('@'):
                    if fileline[1:] in donorreadlist:
                        addedcount = addedcount + 1
                    else:
                        self.vaseutillogger.info("Read " + str(fileline) + " was not added.")
                    next(gzfile)    # Skip the read sequence line
                    next(gzfile)    # Skip the '+' line
                    next(gzfile)    # Skip the read qualities line
        return addedcount
