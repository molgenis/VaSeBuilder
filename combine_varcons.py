import os
import argparse
from VariantContextFile import VariantContextFile


# Returns the command line parameters
def get_parameters():
    combine_varcon_params = argparse.ArgumentParser()
    combine_varcon_params.add_argument("-l", "--listfile", dest="listfile",
                                       help="Listfile containing locations to variant context files to combine")
    combine_varcon_params.add_argument("-i", "--infiles", dest="infiles",
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


# Do the actual work
cli_params = get_parameters()
varcons_to_combine = []
if cli_params['listfile'] is not None:
    varcons_to_combine = read_varcon_list_file(cli_params['listfile'])
else:
    varcons_to_combine = [infile for infile in cli_params['infiles'] if os.path.isfile(infile)]

main_varconfile = VariantContextFile(varcons_to_combine[0])
for varconloc in varcons_to_combine[1:]:
    add_varconfile = VariantContextFile(varconloc)
    main_varconfile.add_variant_context_file(add_varconfile)
main_varconfile.write_variant_context_file(cli_params['outfile'], "")
