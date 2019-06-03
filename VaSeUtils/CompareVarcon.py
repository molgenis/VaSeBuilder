import logging
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class CompareVarcon:
    def __init__(self):
        self.vaseUtilLogger = logging.getLogger("VaSeUtil_Logger")

    # Performs the main analysis
    def main(self, varconfileloc1, varconfileloc2):
        self.vaseUtilLogger.info("Running VaSe util CompareVarcon")
        varconfile1 = VariantContextFile(varconfileloc1)
        varconfile2 = VariantContextFile(varconfileloc2)
        self.compare_varcon_files(varconfile1, varconfile2)
        self.vaseUtilLogger.info("Finished running VaSeUtil CompareVarcon")

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
    def display_shared_varcon_differences(self, varconFile1, varconFile2):
        variantContexts1 = varconFile1.getVariantContextsById()
        variantContexts2 = varconFile2.getVariantContextsById()
        sharedVarconList = varconFile1.get_variant_contexts_intersect(varconFile2)
        sharedVarconDiffs = varconFile1.compare(varconFile2, sharedVarconList)
        diffNumberMap = varconFile1.getVariantContextFields()
        
        # Iterate over the shared variant contexts differences
        print("==================================================")
        print("[-Differing Variant Contexts-]")
        for varconId, varconDiffs in sharedVarconDiffs.items():
            msg= str(varconId)+": "
            for diffield, diffvalues in varconDiffs:
                msg += str(diffNumberMap[diffield])+ " (" +str(diffvalues[0])+ "|" +str(diffvalues[1])+ "), "
            print(msg)
