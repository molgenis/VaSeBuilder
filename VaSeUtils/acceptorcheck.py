import logging
import gzip

"""
Can be used to check that the acceptor reads that should be removed are indeed removed from the fastq file(s).
This is done based on read identifiers.
"""


class AcceptorCheck:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Runs the check to determine if none of the variant context acceptor BAM reads are in the two VaSe produced FastQ
    # files.
    def main(self, readslist, vasefq1, vasefq2):
        self.vaseutillogger.info("Running VaSe util AcceptorCheck")

        # Check if the acceptor reads have indeed not been placed to the VaSe R1 FastQ
        self.vaseutillogger.info("Checking the R1 VaSe FastQ file")
        r1_removed = self.get_number_of_reads_removed(vasefq1, readslist)
        self.vaseutillogger.info("Excluded " + str(r1_removed) + " of " + str(len(readslist))
                                 + " acceptor reads from the R1 VaSe FastQ file.")

        # Check if the acceptor reads have indeed not been added to R2
        self.vaseutillogger.info("Checking the R2 VaSe FastQ file")
        r2_removed = self.get_number_of_reads_removed(vasefq2, readslist)
        self.vaseutillogger.info("Excluded " + str(r2_removed) + " of " + str(len(readslist))
                                 + " acceptor reads from the R2 VaSe FastQ file.")
        
        if r1_removed == len(readslist) and r2_removed == len(readslist):
            self.vaseutillogger.info("All acceptor reads have been excluded from the new VaSe FastQ files")
        else:
            self.vaseutillogger.info("There are still some acceptor reads left in the VaSe FastQ file(s)")
        self.vaseutillogger.info("Finished running VaSe util AcceptorCheck")

    # Returns the number of removed reads for R1/R2 VaSe fastq files.
    def get_number_of_reads_removed(self, vasefqfiles, readidlist):
        removedreads = len(readidlist)
        for vfqFile in vasefqfiles:
            removedreads = self.check_acceptor_reads_removed(vfqFile, readidlist, removedreads)
        return removedreads

    # Checks whether the identified acceptor reads have indeed been removed from a specified VaSe FastQ file.
    def check_acceptor_reads_removed(self, gzresultsfile, acceptorreadlist, removalcount):
        try:
            with gzip.open(gzresultsfile, "rt") as gzfile:
                for fileline in gzfile:
                    fileline = fileline.strip()
                    if fileline.startswith("@"):
                        if fileline[1:] in acceptorreadlist:
                            self.vaseutillogger.info("Read " + str(fileline[1:]) + " was not removed")
                            removalcount = removalcount - 1
                        next(gzfile)    # Skip the sequence line
                        next(gzfile)    # Skip the line with '+'
                        next(gzfile)    # Skip the sequence qualities line
        except IOError:
            self.vaseutillogger.critical(f"Could not open fastq file {gzresultsfile}")
        finally:
            return removalcount
