import logging
import pysam

from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareMappingPositions:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main analysis. Compares the mapping positions of added donor reads to pipeline mapping results
    def main(self, varconfileloc, dlistfileloc, newmappedbamloc):
        donorfiles = self.vuh.read_donorlist_file(dlistfileloc)
        varconfile = VariantContextFile(varconfileloc)
        varcons_bysample = varconfile.get_variant_contexts_by_sampleid()

        try:
            newmappedbamfile = pysam.AlignmentFile(newmappedbamloc)
            for sampleid in varcons_bysample:
                if sampleid in donorfiles:
                    self.process_variant_contexts(varcons_bysample[sampleid], newmappedbamfile, donorfiles[sampleid])
                else:
                    self.vaseutillogger.warning(f"No valid BAM file found for sample {sampleid}")
        except IOError:
            self.vaseutillogger.warning(f"Could not open BAM file: {newmappedbamloc}")

    # Processes a variant context file.
    def process_variant_contexts(self, variantcontexts, newmappedbamfile, bamfileloc):
        try:
            bamfile = pysam.AlignmentFile(bamfileloc)
            for varcon in variantcontexts:
                search_chrom = varcon.get_variant_context_chrom()
                search_start = varcon.get_variant_context_start()
                search_end = varcon.get_variant_context_end()
                readidfilter = varcon.get_donor_reads()

                donordata  = self.get_bam_mapping_positions(bamfile, search_chrom, search_start, search_start,
                                                            readidfilter)
                mappeddata = self.get_bam_mapping_positions(newmappedbamfile, search_chrom, search_start, search_end,
                                                            readidfilter)

                self.display_differences(donordata, mappeddata)
        except IOError:
            self.vaseutillogger.warning(f"Could not open BAM file: {bamfileloc}")

    # Returns the left positions and CIGAR strings of reads.
    # readidfilter should contain the donor read ids of a variant context
    def get_bam_mapping_positions(self, bamfile, searchchrom, searchstart, searchend, readidfilter):
        mapping_positions = {}
        for bamread in bamfile.fetch(searchchrom, searchstart, searchend):
            if bamread.query_name in readidfilter:
                mapping_positions[bamread.query_name] = (bamread.reference_start, bamread.cigarstring)
        return mapping_positions

    def display_differences(self,donorpositions, newmappositions):
