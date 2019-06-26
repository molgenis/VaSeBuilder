import logging
import pysam
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile
from VcfVariant import VcfVariant

"""Checks whether VCF variants that were not used to construct a variant context are located in one"""


class VcfVariantsInVariantContexts:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main analysis
    def main(self, vcf_file, varcon_file):
        vcfvars_incontexts = self.check_nonvarcons_variant_in_contexts(vcf_file, varcon_file)
        if len(vcfvars_incontexts['is_context']) > 0 or len(vcfvars_incontexts['in_context']) > 0:
            self.display_results(vcfvars_incontexts)
        else:
            print(f"VCF file {vcf_file} has no variants that are a context or are located in a context")

    # Checks the non variantcontext VCF variants
    def check_nonvarcons_variant_in_contexts(self, vcf_fileloc, varconfile):
        vcfvariants_incontexts = {'is_context': [], 'in_context': []}
        try:
            vcf_file = pysam.VariantFile(vcf_fileloc, "r")
            for vcfvariant in vcf_file.fetch():
                varianttype = self.get_variant_type(vcfvariant.ref, vcfvariant.alts)
                variant_startstop = self.determine_indel_range(vcfvariant.pos, vcfvariant.ref, vcfvariant.alts)
                vcfvardata = VcfVariant(vcfvariant.chrom, vcfvariant.pos, vcfvariant.ref, vcfvariant.alts,
                                        vcfvariant.filter, varianttype)

                # Check if the VCF variant is a context or within a context
                if varconfile.has_variant_context(f"{vcfvariant.chrom}_{vcfvariant.pos}"):
                    vcfvariants_incontexts['is_context'].append(vcfvardata)
                elif varconfile.variant_is_in_context(varianttype, vcfvariant.chrom, variant_startstop[0],
                                                      variant_startstop[1]):
                    vcfvariants_incontexts['in_context'].append(vcfvardata)
            vcf_file.close()
        except IOError:
            self.vaseutillogger.warning(f"Could not open VCF file {vcf_fileloc}")
        finally:
            return vcfvariants_incontexts

    # Displays the results from the check_nonvarcons function
    def display_results(self, context_results):
        if len(context_results['is_context']) > 0:
            print("=====VCF VARIANTS THAT ARE A CONTEXT=====")
            self.display_variants(context_results['is_context'])
            print("")
        if len(context_results['in_context']) > 0:
            print("=====VCF VARIANTS THAT ARE IN A CONTEXT=====")
            self.display_variants(context_results['in_context'])

    # Returns the type of the variant needed for checking whether the variant is in a context
    def get_variant_type(self, vcfvar_ref, vcfvar_alts):
        maxreflength = max([len(x) for x in vcfvar_ref.split(",")])
        maxaltlength = max([len(x) for x in vcfvar_alts])
        if maxreflength > 1 and maxaltlength > 1:
            return "indel"
        return "snp"

    # Determines the range of the indel
    def determine_indel_range(self, variantpos, variantref, variantalts):
        indelstart = variantpos
        indelstop = variantpos + max(max([len(x) for x in variantref.split(",")]), max([len(x) for x in variantalts]))
        return [indelstart, indelstop]

    # Print the variants in a variantlist
    def display_variants(self, variantlist):
        for vcfvar in variantlist:
            print(f"{vcfvar.to_string()}")
