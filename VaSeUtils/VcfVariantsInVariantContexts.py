import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

"""Checks whether VCF variants that were not used to construct a variant context are located in one"""


class VcfVariantsInVariantContexts:
    def __init__(self, vaseutilhelper):
        self.vaseutillogger = logging.getLogger("VaSe_Logger")
        self.vuh = vaseutilhelper

    # Performs the main analysis
    def main(self, varconfile, vcffilelist):
        print("aap")

    # Checks the non variantcontext VCF variants
    def check_nonvarcon_variant_in_context(self):
        print("aap")
