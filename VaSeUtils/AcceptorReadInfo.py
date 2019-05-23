import logging
import pysam

# Import required VaSeUtils classes
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class AcceptorReadInfo:
    def __init__(self, vaseuhelper):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.vuh = vaseuhelper
    
    
    # Performs all the analysis steps
    def main(self, acceptorBamFile, vcFileLoc, sampleFilter=None, varconFilter=None, readIdFilter=None):
        self.vaseutillogger.info("Running VaSe util AcceptorReadInfo")
        varconfile = VariantContextFile(vcFileLoc, sampleFilter, varconFilter)
        self.get_acceptor_read_info(acceptorBamFile, varconfile)
        self.vaseutillogger.info("Finished running VaSe util AcceptorReadInfo")

    # Obtains the read info for the selected variant context acceptor reads (all reads if no filters were set)
    def get_acceptor_read_info(self, acceptorbam, varconfile, readidfilter=None):
        try:
            abamfile = pysam.AlignmentFile(acceptorbam)
            for variantcontext in varconfile.get_variant_contexts():
                search_chrom = variantcontext.get_variant_context_chrom()
                search_start = variantcontext.get_variant_context_start()
                search_stop = variantcontext.get_variant_context_end()

                if search_chrom and search_start and search_stop:
                    print(f"Variant context {variantcontext.get_variant_context_id()} ;; from"
                          f"sample {variantcontext.get_variant_context_sample()}")
                    self.get_read_info(variantcontext.get_acceptor_reads(), search_chrom, search_start, search_stop,
                                       abamfile, readidfilter)
        except IOError:
            self.vaseutillogger.critical(f"Could not open acceptor BAM file: {acceptorbam}")

    # Obtain info for a specific list of variant context acceptor reads
    def get_read_info(self, readslist, searchchrom, searchstart, searchend, bamfile, readidfilter=None):
        for bamread in bamfile.fetch(searchchrom, searchstart, searchend):
            if bamread.query_name in readslist and self.vuh.passes_filter(bamread.query_name, readidfilter):
                print(bamread.to_string())
                print(self.get_mate_info(bamfile, bamread))

    # Obtain the read mate info
    def get_mate_info(self, bamfile, bamread):
        try:
            materead = bamfile.mate(bamread)
            return materead.to_string()
        except ValueError:
            return ""
