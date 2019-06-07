import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareVarcon:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")

    # Performs the main analysis
    def main(self, varconfileloc1, varconfileloc2):
        self.vaseutillogger.info("Running VaSe util CompareVarcon")
        varconfile1 = VariantContextFile(varconfileloc1)
        varconfile2 = VariantContextFile(varconfileloc2)
        self.compare_varcon_files(varconfile1, varconfile2)
        self.vaseutillogger.info("Finished running VaSeUtil CompareVarcon")

    # Performs the comparison between two variant context files.
    def compare_varcon_files(self, varconfile1, varconfile2):
        self.display_general_info(varconfile1, varconfile2)
        self.display_shared_varcon_differences(varconfile1, varconfile2)

    # Displays general info about the two variant contexts files
    def display_general_info(self, varconfile1, varconfile2):
        print("==================================================")
        print("[-Variant Context General Info-]")
        print("\tVarcon1\tVarcon2")
        print(f"#Contexts:\t{varconfile1.get_number_of_contexts()}\t{varconfile2.get_number_of_contexts()}")
        print(f"#Difference:\t{varconfile1.get_variant_contexts_difference(varconfile2)}"
              f"\t{varconfile2.get_variant_contexts_difference(varconfile1)}")
        print("--------------------")
        print(f"#Union: {varconfile1.get_variant_contexts_union(varconfile2)}")
        print(f"#Intersect: {varconfile1.get_variant_contexts_intersect(varconfile2)}")
        print(f"#Symmetric Difference: {varconfile1.get_variant_contexts_symmetric_difference(varconfile2)}")

    # Displays which variant contexts are in both but differ (context start/end, etc)
    def display_shared_varcon_differences(self, varconfile1, varconfile2):
        #variant_contexts1 = varconfile1.getVariantContextsById()
        #variant_contexts2 = varconfile2.getVariantContextsById()
        shared_varcon_list = varconfile1.get_variant_contexts_intersect(varconfile2)
        shared_varcon_diffs = varconfile1.compare(varconfile2, shared_varcon_list)
        diff_number_map = varconfile1.get_variant_context_fields()
        
        # Iterate over the shared variant contexts differences
        print("==================================================")
        print("[-Differing Variant Contexts-]")
        for varconid, varcondiffs in shared_varcon_diffs.items():
            msg = f"{varconid}: "
            for diffield, diffvalues in varcondiffs:
                msg += f"{diff_number_map[diffield]} ({diffvalues[0]}|{diffvalues[1]}), "
            print(msg)
