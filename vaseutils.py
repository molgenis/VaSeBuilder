import argparse

from VariantContextFile import VariantContextFile
from VaSeUtils.UtilParamCheck import UtilParamCheck
from VaSeUtils.VaSeUtilRunner import VaSeUtilRunner


def get_util_parameters(utilparamcheck):
    """Gathers and returns the set parameters

    Returns
    -------
    dict
        Parameter names and the set values
    """
    utils_argparse = argparse.ArgumentParser()
    utils_argparse.add_argument("-u", "--util", dest="util", required=True,
                                choices=utilparamcheck.get_valid_vaseutils(), help="Util to run")
    utils_argparse.add_argument("-vcf", "--vcflist", dest="vcflist", help="VCF list file")
    utils_argparse.add_argument("-vci", "--varconin", dest="varconin", help="Variant context input file")
    utils_argparse.add_argument("-vl", "--variantlist", dest="variantlist", help="Variant list file")
    utils_argparse.add_argument("-df", "--donorfiles", dest="donorfiles", help="Donor list file")
    return vars(utils_argparse.parse_args())


def run_selected_util(util_to_run, util_parameters, vase_util_runner, vase_util_paramcheck):
    """Runs the selected VaSe util with the provided parameters.

    Parameters
    ----------
    util_to_run : str
        VaSeUtil to run
    util_parameters : dict
        Parameters and values to run selected util with
    vase_util_runner : VaSeUtilRunner
        Object that runs a selected VaSeUtil
    vase_util_paramcheck : UtilParamCheck
        Utility to check parameters
    """
    # Check whether to read a provided VariantContextFile
    varconfile = None
    if util_parameters["varconin"] is not None:
        varconfile = VariantContextFile(util_parameters["varconin"])

    if util_to_run == "checkdonorfiles":
        vase_util_runner.check_donor_files_exist(util_parameters["donorfiles"])
    if util_to_run == "loginfo":
        if vase_util_paramcheck.optional_parameter_is_set("logfilter", util_parameters):
            vase_util_runner.log_info(util_parameters["vaselog"], util_parameters["logfilter"])
        else:
            vase_util_runner.log_info(util_parameters["vaselog"])
    if util_to_run == "subvcfvarcon":
        vase_util_runner.subset_vcf_by_variant_contexts(varconfile, util_parameters["vcflist"])
    if util_to_run == "subvcfvbarlist":
        vase_util_runner.subset_vcf_by_variant_list(util_parameters["variantlist"], util_parameters["vcflist"])


# Run a selected VaSeUtil
if __name__ == "__main__":
    vur = VaSeUtilRunner()
    upc = UtilParamCheck()

    utilparams = get_util_parameters(upc)
    if utilparams["util"] is not None:
        if upc.required_params_set(utilparams["util"], utilparams):
            run_selected_util(utilparams["util"], utilparams, vur, upc)
