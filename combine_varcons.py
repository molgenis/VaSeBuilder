import os
import argparse
from VariantContextFile import VariantContextFile


# Returns the command line parameters
def get_parameters():
    combine_varcon_params = argparse.ArgumentParser()
    combine_varcon_params.add_argument("-l", "--listfile", dest="listfile",
                                       help="Listfile containing locations to variant context files to combine")
    combine_varcon_params.add_argument("-i", "--infiles", dest="infiles", nargs="+",
                                       help="Variant context files to combine.")
    combine_varcon_params.add_argument("-o", "--out", dest="outfile")
    combine_varcon_args = vars(combine_varcon_params.parse_args())
    return combine_varcon_args


# Reads the variant context list file
def read_varcon_list_file(varcon_listfile):
    varcon_files = []
    try:
        with open(varcon_listfile, 'r') as vclf:
            for fileline in vclf:
                if os.path.isfile(fileline.strip()):
                    varcon_files.append(fileline.strip())
    except IOError:
        print(f"Could not read variant context list file: {varcon_listfile}")
    finally:
        return varcon_files


def combine_varcons(varconfiles):
    varcondata = {}
    for varconfileloc in varconfiles:
        varcondata = process_varconfile(varconfileloc, varcondata)
    return varcondata


def process_varconfile(varconfileloc, varcondata):
    with open(varconfileloc, 'r') as vrcfile:
        for fileline in vrcfile:
            if not fileline.startswith("#"):
                varconline = fileline.strip()
                filelinedata = fileline.strip().split("\t")

                if not detect_collision(filelinedata, varcondata):
                    varcondata[filelinedata[0]] = varconline
        return varcondata


def detect_collision(varconentry, varcondatamap):
    if varconentry[0] in varcondatamap:
        return True
    else:
        for varcon in varcondatamap.values():
            if varcon[2] == varconentry[2]:
                if int(varcon[4]) <= int(varconentry[5]) and int(varconentry[4]) <= int(varcon[5]):
                    return True
        return False


def write_new_variantcontextfile(outfileloc, varcondatamap):
    with open(outfileloc, 'w') as vrcoutfile:
        vrcoutfile.write("#ContextId\tDonorSample\tChrom\tOrigin\tStart\tEnd\tAcceptorContextLength\t"
                         "DonorContextLength\tAcceptorReads\tDonorReads\tADratio\tAcceptorReadsIds\tDonorReadIds\n")
        for varcondataentry in varcondatamap.values():
            vrcoutfile.write(f"{varcondataentry}\n")


# Do the actual work
cli_params = get_parameters()
varcons_to_combine = []
if cli_params['listfile'] is not None:
    varcons_to_combine = read_varcon_list_file(cli_params['listfile'])
else:
    varcons_to_combine = [infile for infile in cli_params['infiles'] if os.path.isfile(infile)]

varcondata = combine_varcons(varcons_to_combine)
write_new_variantcontextfile(cli_params['outfile'], varcondata)

"""
main_varconfile = VariantContextFile(varcons_to_combine[0])
for varconloc in varcons_to_combine[1:]:
    add_varconfile = VariantContextFile(varconloc)
    main_varconfile.add_variant_context_file(add_varconfile)
main_varconfile.write_variant_context_file(cli_params['outfile'], "")
"""
