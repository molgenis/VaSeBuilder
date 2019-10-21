import argparse

from VariantContextFile import VariantContextFile
from VaSeUtils.UtilParamCheck import UtilParamCheck
from VaSeUtils.VaSeUtilHelper import VaSeUtilHelper
from VaSeUtils.VaSeUtilRunner import VaSeUtilRunner


def get_util_parameters():
    """Gathers and returns the set parameters

    Returns
    -------
    dict
        Parameter names and the set values
    """
    utils_argparse = argparse.ArgumentParser()
    utils_argparse.add_argument("-u", "--util", dest="util", help="")
    utils_argparse.add_argument("-vcf", "--vcflist", dest="vcflist", help="")
    utils_argparse.add_argument("-vci", "--varconin", dest="varconin", help="")
    utils_argparse.add_argument("-vl", "--variantlist", dest="variantlist", help="")
    utils_argparse.add_argument("-df", "--donorfiles", dest="donorfiles", help="")
    return vars(utils_argparse)


def run_selected_util(util_to_run, util_parameters):
    """Runs the selected VaSe util with the provided parameters.

    Parameters
    ----------
    util_to_run : str
        VaSeUtil to run
    util_parameters : dict
        Parameters and values to run selected util with
    """
    # Check whether to read a provided VariantContextFile
    varconfile = None
    if util_parameters["varconin"] is not None:
        varconfile = VariantContextFile(util_parameters["varconin"])

    if util_to_run == "subvcfvarcon":
        vur.subset_vcf_by_variant_contexts(varconfile, util_parameters["vcflist"])
    if util_to_run == "subvcfvbarlist":
        vur.subset_vcf_by_variant_list(util_parameters["variantlist"], util_parameters["vcflist"])
    if util_to_run == "checkdonorfiles":
        vur.check_donor_files_exist(util_parameters["donorfiles"])


# Run a selected VaSeUtil
utilparams = get_util_parameters()

vur = VaSeUtilRunner()
upc = UtilParamCheck()
if utilparams["util"] is not None:
    if upc.required_params_set(utilparams["util"], utilparams):
        run_selected_util(utilparams["util"], utilparams)
