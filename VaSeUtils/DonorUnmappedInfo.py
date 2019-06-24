import logging
import pysam
from DonorBamRead import DonorBamRead


class DonorUnmappedInfo:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main functionality
    def main(self, unmappedfileloc, variantcontexts, donorbamlocs):
        unmapped_data = self.read_donor_unmapped_matefile(unmappedfileloc)
        for sampleid in unmapped_data:

    # Reads the donor bam unmapped file
    def read_donor_unmapped_matefile(self, unmappedfileloc):
        donor_unmapped_ids = {}
        try:
            with open(unmappedfileloc, 'r') as unmappedfile:
                next(unmappedfile)    # Skip the header line
                for fileline in unmappedfile:
                    filelinedata = fileline.strip().split("\t")
                    if filelinedata[1] not in donor_unmapped_ids:
                        donor_unmapped_ids[filelinedata[1]] = {}
                    donor_unmapped_ids[filelinedata[1]][filelinedata[0]] = filelinedata[2].split(";")
        except IOError:
            self.vaseutillogger.critical(f"Could not open unmapped domor mate file: {unmappedfileloc}")
        finally:
            return donor_unmapped_ids

    # Fetches the unmapped mate info for one sample from the donor BAM file
    def fetch_sample_unmapped_info(self, unmapped_sample_data, donorbamloc):
        sample_unmapped_info = {}
        try:
            dbamfile = pysam.AlignmentFile(donorbamloc, "rb")
            dbamfile.close()
        except IOError:
            self.vaseutillogger.critical(f"Could not open donor BAM file {donorbamloc}")
        finally:
            return sample_unmapped_info

    # Fetches the unmapped mate info for one context from a donor BAM file
    def fetch_context_unmapped_info(self, unmapped_context_data, donorbamfile):
        for fetched_read in donorbamfile.fetch()
