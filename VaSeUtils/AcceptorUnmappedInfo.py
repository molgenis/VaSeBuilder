import logging
import pysam
from DonorBamRead import DonorBamRead


class AcceptorUnmappedInfo:
    # Constructor that sets the logger and VaSeUtilHelper
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main analysis
    def main(self, unmappedfileloc, variantcontexts, acceptorbamloc):
        unmapped_mates = self.vuh.read_unmapped_matefile(unmappedfileloc)
        try:
            abamfile = pysam.AlignmentFile(acceptorbamloc)
            for contextid in variantcontexts:
                search_chrom = variantcontexts[contextid].get_variant_context_chrom()
                search_start = variantcontexts[contextid].get_variant_context_start()
                search_end = variantcontexts[contextid].get_variant_context_end()

                unmapped_readdata = self.get_unmmapped_mate_info(search_chrom, search_start, search_end, abamfile,
                                                                 unmapped_mates[contextid])
                self.display_unmapped_read_info(unmapped_readdata)
            abamfile.close()
        except IOError:
            self.vaseutillogger.critical(f"Could not open acceptor BAM file {acceptorbamloc}")
        finally:
            return unmapped_mates

    # Fetches unmapped read info fro the BAM file
    def get_unmmapped_mate_info(self, searchchrom, searchstart, searchend, abamfile, readidfilter):
        fetch_data = {'mapped': [], 'unmapped': []}
        for fetch_read in abamfile.fetch(searchchrom, searchstart, searchend):
            if self.vuh.passes_filter(fetch_read.query_name, readidfilter):
                dbr_read = DonorBamRead(fetch_read.query_name, self.vuh.get_read_pair_num(), fetch_read.reference_name,
                                        fetch_read.reference_start, fetch_read.infer_read_length(),
                                        fetch_read.get_forward_sequence(),
                                        "".join([chr((x + 33)) for x in fetch_read.get_forward_qualities()]),
                                        fetch_read.mapping_quality)

                # Check whether to add to the mapped or unmapped list
                if fetch_read.is_unmapped:
                    fetch_data['unmapped'].append(dbr_read)
                else:
                    fetch_data['mapped'].append(dbr_read)
        return fetch_data

    # Displays the unmapped read info for a context
    def display_unmapped_read_info(self, unmapped_readdata):
        mapped_reads = sorted(unmapped_readdata['mapped'], key=lambda x: x.get_bam_read_id(), reverse=False)
        unmapped_reads = sorted(unmapped_readdata['unmapped'], key=lambda x: x.get_bam_read_id(), reverse=False)

        # Display the read pairs (mapped first, unmapped second)
        for x in range(0, len(mapped_reads)):
            print(f"{mapped_reads[x].to_string()}")
            print(f"{unmapped_reads[x].to_string()}")
