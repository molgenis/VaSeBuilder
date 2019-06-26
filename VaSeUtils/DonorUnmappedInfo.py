import logging
import pysam
from DonorBamRead import DonorBamRead


class DonorUnmappedInfo:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main functionality
    def main(self, unmappedfileloc, variantcontexts, donorbamlocs):
        donorbams = self.vuh.read_donor_list_file(donorbamlocs)
        unmapped_data = self.read_donor_unmapped_matefile(unmappedfileloc)
        for sampleid in unmapped_data:
            if sampleid in donorbams:
                contextids = list(unmapped_data[sampleid].keys())
                self.fetch_sample_unmapped_info(unmapped_data, variantcontexts, contextids, donorbams[sampleid])

    # Reads the donor bam unmapped file
    def read_donor_unmapped_matefile(self, unmappedfileloc):
        donor_unmapped_ids = {}
        try:
            with open(unmappedfileloc, "r") as unmappedfile:
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
    def fetch_sample_unmapped_info(self, unmapped_sample_data, variantcontexts, contextidlist, donorbamloc):
        sample_unmapped_info = {}
        try:
            dbamfile = pysam.AlignmentFile(donorbamloc, "rb")
            for contextid in contextidlist:
                searchchrom = variantcontexts[contextid].get_variant_context_chrom()
                searchstart = variantcontexts[contextid].get_variant_context_start()
                searchend = variantcontexts[contextid].get_variant_context_end()
                context_umapdata = self.fetch_context_unmapped_info(searchchrom, searchstart, searchend, dbamfile,
                                                                    unmapped_sample_data[contextid])
                self.display_unmapped_read_info(context_umapdata)
            dbamfile.close()
        except IOError:
            self.vaseutillogger.critical(f"Could not open donor BAM file {donorbamloc}")
        finally:
            return sample_unmapped_info

    # Fetches the unmapped mate info for one context from a donor BAM file
    def fetch_context_unmapped_info(self, searchchrom, searchstart, searchend, donorbamfile, readidfilter):
        fetch_data = {'mapped': [], 'unmapped': []}
        for fetch_read in donorbamfile.fetch(searchchrom, searchstart, searchend):
            if self.vuh.passes_filter(fetch_read.query_name, readidfilter):
                dbr_read = DonorBamRead(fetch_read.query_name, self.vuh.get_read_pair_num(), fetch_read.reference_name,
                                        fetch_read.reference_start, fetch_read.infer_read_length(),
                                        fetch_read.get_forward_sequence(),
                                        "".join([chr((x + 33)) for x in fetch_read.get_forward_qualities()]),
                                        fetch_read.mapping_quality)
                if fetch_read.is_unmapped:
                    fetch_data['unmapped'].append(dbr_read)
                else:
                    fetch_data['mapped'].append(dbr_read)
        return fetch_data

    # Display the actual data of the mapped and unmapped reads per pair
    def display_unmapped_read_info(self, readdata):
        mapped_reads = sorted(readdata['mapped'], key=lambda x: x.get_bam_read_id(), reverse=False)
        unmapped_reads = sorted(readdata['unmapped'], key=lambda x: x.get_bam_read_id(), reverse=False)

        for y in range(0, len(mapped_reads)):
            print(f"{mapped_reads[y].to_string()}")
            print(f"{unmapped_reads[y].to_string()}")
