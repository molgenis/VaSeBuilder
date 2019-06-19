import logging
import gzip
from donorcheck import DonorCheck
from acceptorcheck import AcceptorCheck


class CheckVaSeFastQ:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Checks that the VaSe R1 and R2 FastQ have the correct number of lines.
    def main(self, templatefq1, vasefq1, templatefq2, vasefq2, acceptorreadlist, donorreadlist):
        self.vaseutillogger.info("Running VaSe util CheckVaSeFastQ")
        self.check_vase_fastq(templatefq1, vasefq1, templatefq2, vasefq2, acceptorreadlist, donorreadlist)
        self.vaseutillogger.info("Finished running VaSe util CheckVaSeFastQ")

    # Checks the VaSe produced FastQ files (is the size as should be, are all required template reads added to the VaSe
    # produced fastq files)
    def check_vase_fastq(self, templatefq1, vasefq1, templatefq2, vasefq2, acceptorreadlist, donorreadlist):
        # Check the VaSe produced R1 fastq file.
        if self.check_vase_fastq_size(templatefq1, vasefq1, acceptorreadlist, donorreadlist):
            self.vaseutillogger.info(f"VaSe produced FastQ R1 file {vasefq1} seems to have the expected size.")
        else:
            self.vaseutillogger.info(f"VaSe produced FastQ R1 file {vasefq1} does not seem to have the expected size.")
        
        # Check the VaSe produced R2 fastq file.
        if self.check_vase_fastq_size(templatefq2, vasefq2, acceptorreadlist, donorreadlist):
            self.vaseutillogger.info(f"VaSe produced FastQ R2 file {vasefq2} seems to have the expected size.")
        else:
            self.vaseutillogger.info(f"VaSe produced FastQ R2 file {vasefq2} does not seem to have the expected size.")
        
        # Lastly check that the reads of the VaSe FastQ and acceptor reads are the same as those in the template
        # (original) fastq files
        #checkReadsInFastQs()    # Check for R1
        #checkReadsInFastQs()    # Check for R2

    # Compares an original (template) and VaSe produced fastq file
    def check_vase_fastq_size(self, templatefq, vasefq, donorreads, acceptorreads):
        acheck = AcceptorCheck()
        dcheck = DonorCheck()
        numof_acceptorreads = acheck.check_acceptor_reads_removed(vasefq, acceptorreads)
        numof_donorreads = dcheck.check_donor_reads_added(vasefq, donorreads)
        
        # Get the number of lines for the template/original fastq file
        tmplfq_num_of_lines = 0
        if templatefq.endswith(".fastq") or templatefq.endswith(".fq"):
            tmplfq_num_of_lines = self.get_fastq_length(templatefq)
        else:
            tmplfq_num_of_lines = self.get_fastq_gz_length(templatefq)
        
        # Get the number of lines for the vase created fastq file
        vasefq_num_of_lines = 0
        if vasefq.endswith(".fastq") or vasefq.endswith(".fq"):
            vasefq_num_of_lines = self.get_fastq_length(vasefq)
        else:
            vasefq_num_of_lines = self.get_fastq_gz_length(vasefq)
        
        # Calculate what the length of the VaSe FastQ should be, multiply with 4 as each read has four lines.
        control_vase_length = ((tmplfq_num_of_lines - (numof_acceptorreads*4)) + (numof_donorreads*4))
        
        # Check whether the the calculated length (control_vase_length) and the actual length of the vase fastq are the
        # same. Also check calculated fq length + 1 in case the actual files have an empty line at the end of the file.
        if control_vase_length == vasefq_num_of_lines or (control_vase_length+1) == vasefq_num_of_lines:
            return True
        return False

    # Returns the number of lines in a plain text fastq file.
    # Thanks to SilentGhost at https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    def get_fastq_length(self, fqfileloc):
        linenum = -1
        with open(fqfileloc, "r") as fqfile:
            for linenum, fileline in enumerate(fqfile):
                pass
        return linenum + 1

    # Returns the number of lines in a gzipped fastq file.
    def get_fastq_gz_length(self, fqfileloc):
        linenum = -1
        with gzip.open(fqfileloc, "rt") as fqfile:
            for linenum, fileline in enumerate(fqfile):
                pass
        return linenum + 1
