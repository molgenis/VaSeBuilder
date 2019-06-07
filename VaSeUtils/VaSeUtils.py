import argparse
import logging

# Import the VaSeUtils classes
from acceptorcheck import AcceptorCheck
from AcceptorReadInfo import AcceptorReadInfo
from checkVaseFastQ import CheckVaSeFastQ
from CompareAcceptor import CompareAcceptor
from CompareAcceptorContexts import CompareAcceptorContexts
from CompareDonorContexts import CompareDonorContexts
from CompareVarcon import CompareVarcon
from donorcheck import DonorCheck
from DonorReadInfo import DonorReadInfo
from LogInfo import LogInfo
from UtilParamCheck import UtilParamCheck
from VaSeUtilHelper import VaSeUtilHelper
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class VaSeUtils:
    def __init__(self):
        self.upc = UtilParamCheck()
        self.vuh = VaSeUtilHelper()
        self.valid_utils = ['acceptorcheck', 'acceptorreadinfo', 'checkdonorfiles', 'checkfastq', 'compareacceptor',
                            'comparedonor', 'comparefastq', 'comparevarcon', 'donorcheck', 'donorreadinfo', 'loginfo',
                            'subsetvarcon', 'unmappedinfo', 'varcondata', 'vasecompare']
        self.vaseutillogger = self.start_logger()

    # Runs all specified VaSe utils
    def main(self):
        vaseu_args = self.get_vaseutils_parameters()
        self.vaseutillogger.info("Running selected VaSeUtil program(s)")
        
        # Run one or more selected utils
        for utilToRun in vaseu_args['util']:
            self.run_selected_util(utilToRun, vaseu_args)
        self.vaseutillogger.info("Ran selected VaSeUtil program(s)")

    # Returns the set parameter values.
    def get_vaseutils_parameters(self):
        vaseutil_argpars = argparse.ArgumentParser(description="Run a specific VaSe Util program")
        vaseutil_argpars.add_argument("-u", "--util", nargs="+", dest="util", choices=self.valid_utils, required=True,
                                      help="The utility to run.")
        vaseutil_argpars.add_argument("-l", "--log", dest='log', help="Location to write VaSeUtils log file to")
        vaseutil_argpars.add_argument("-df", "--donorfiles", dest="donorfiles", help="File containing the list of used"
                                                                                     "donor VCF/BAM files")
        vaseutil_argpars.add_argument("-vf1", "--vasefq1", dest="vasefq1", nargs="*", help="The VaSe produced R1 FastQ"
                                                                                           "file(s)")
        vaseutil_argpars.add_argument("-vf2", "--vasefq2", dest="vasefq2", nargs="*", help="The VaSe produced R2 FastQ"
                                                                                           "file(s)")
        vaseutil_argpars.add_argument("-ov1", "--othervasefq1", dest="othervasefq1", nargs="*", help="The other VaSe"
                                                                                                     "produced R1 FastQ"
                                                                                                     "file(s)")
        vaseutil_argpars.add_argument("-ov2", "--othervasefq2", dest="othervasefq2", nargs="*", help="The other VaSe"
                                                                                                     "produced R2 FastQ"
                                                                                                     "file(s)")
        vaseutil_argpars.add_argument("-vl", "--vaselog", dest="vaselog", help="Location to the log file produced by"
                                                                               "VaSeBuilder")
        vaseutil_argpars.add_argument("-tf1", "--templatefq1", dest="templatefq1", nargs="*", help="Template R1 FastQ"
                                                                                                   "file used to"
                                                                                                   "produce the VaSe R1"
                                                                                                   "FastQ file")
        vaseutil_argpars.add_argument("-tf2", "--templatefq2", dest="templatefq2", nargs="*", help="Template R2 FastQ"
                                                                                                   "file used to"
                                                                                                   "produce the VaSe R2"
                                                                                                   "FastQ file")
        vaseutil_argpars.add_argument("-ab", "--acceptorbam", dest="acceptorbam", help="BAM file used as acceptor")
        vaseutil_argpars.add_argument("-vc", "--varcon", dest="varcon", help="VaSe produced variant context file")
        vaseutil_argpars.add_argument("-vc2", "--varcon2", dest="varcon2", help="Other VaSe produced variant context"
                                                                                "file")
        vaseutil_argpars.add_argument("-um", "--unmappedmates", dest="unmappedmates", help="VaSe produced file with"
                                                                                           "read identifers that have"
                                                                                           "unmapped mates")
        vaseutil_argpars.add_argument("-um2", "--unmappedmates2", dest="unmappedmates2", help="Other VaSe produced file"
                                                                                              "with read identifiers"
                                                                                              "that have unmapped"
                                                                                              "mates")
        vaseutil_argpars.add_argument("-sf", "--samplefilter", dest="samplefilter", help="List of sample identifiers to"
                                                                                         "include. Will use all samples"
                                                                                         "if not set")
        vaseutil_argpars.add_argument("-cf", "--chromfilter", dest="chromfilter", help="List of chromosomes to use."
                                                                                       "Will use all chromosomes if not"
                                                                                       "set")
        vaseutil_argpars.add_argument("-pf", "--posfilter", dest="posfilter", help="List of start-end position ranges"
                                                                                   "to use. Will use all positions if"
                                                                                   "not set")
        vaseutil_argpars.add_argument("-vf", "--varconfilter", dest="varconfilter", help="List of variant context to"
                                                                                         "use. Will use all variant"
                                                                                         "contexts if not set")
        vaseutil_argpars.add_argument("-lf", "--logfilter", dest="logfilter", help="Filter for which log fields to show"
                                                                                   "(e.g. INFO, DEBUG, WARNING)")
        vaseutil_argpars.add_argument("-rif", "--readidfilter", dest="readidfilter", help="Filter for which reads to"
                                                                                          "obtain info for")
        vaseutil_argpars.add_argument("-i1", "--infile1", dest="infile1", help="First input file")
        vaseutil_argpars.add_argument("-i2", "--infile2", dest="infile2", help="Second input file")
        vaseutil_argpars.add_argument("-o", "--out", dest="outfile", help="Location of the output file")
        return vars(vaseutil_argpars.parse_args())

    # Runs a selected util
    def run_selected_util(self, utiltorun, programparams):
        if self.upc.required_params_set(utiltorun, programparams):
            # Run the AcceptorCheck util.
            if utiltorun == "acceptorcheck":
                varcon_file = VariantContextFile(programparams['varcon'])
                bamread_list = varcon_file.get_all_variant_context_acceptor_read_ids()
                acheck = AcceptorCheck()
                acheck.main(bamread_list, programparams['vasefq1'], programparams['vasefq2'])
            
            # Run the AcceptorReadInfo util.
            if utiltorun == "acceptorreadinfo":
                ari = AcceptorReadInfo(self.vuh)
                ari.main(programparams['acceptorbam'], programparams['varcon'], programparams['samplefilter'],
                         programparams['varconfilter'], programparams['readidfilter'])
            
            # Run the CheckFastQ util.
            if utiltorun == "checkfastq":
                varcon_file = VariantContextFile(programparams['varcon'])
                acceptorreadlist = varcon_file.get_all_variant_context_acceptor_read_ids()
                donorreadlist = varcon_file.get_all_variant_context_donor_read_ids()
                checkf = CheckVaSeFastq()
                checkf.main(programparams['templatefq1'], programparams['vasefq1'], programparams['templatefq2'],
                            programparams['vasefq2'], acceptorreadlist, donorreadlist)

            # Run the CheckDonorFiles util
            if utiltorun == "checkdonorfiles":
                check_dfiles = CheckDonorFilesExist(self.vuh)
                check_dfiles.main(programparams['donorfiles'], programparams['samplefilter'])

            # Run the CompareAcceptorContexts util
            if utiltorun == "compareacceptor":
                cacs = CompareAcceptorContexts(self.vuh)
                cacs.main(programparams['infile1'], programparams['infile2'], programparams['samplefilter'],
                          programparams['varconfilter'], programparams['chromfilter'])

            # Run the CompareDonorContexts util
            if utiltorun == "comparedonor":
                cdcs = CompareDonorContexts(self.vuh)
                cdcs.main(programparams['infile1'], programparams['infile2'], programparams['samplefilter'],
                          programparams['varconfilter'], programparams['chromfilter'])

            # Run the CompareVariantContexts util
            if utiltorun == "comparevarcon":
                cvcs = CompareVarcon()
                cvcs.main(programparams['varcon'], programparams['varcon2'], programparams['samplefilter'],
                          programparams['varconfilter'], programparams['chromfilter'])

            # Run the DonorCheck util.
            if utiltorun == "donorcheck":
                varcon_file = VariantContextFile(programparams['varcon'])
                bamread_list = varcon_file.get_all_variant_context_donor_read_ids()
                dcheck = DonorCheck()
                dcheck.main(bamread_list, programparams['vasefq1'], programparams['vasefq2'])
            
            # Run the DonorReadInfo util.
            if utiltorun == "donorreadinfo":
                dri = DonorReadInfo(self.vuh)
                dri.main(programparams['donorfiles'], programparams['varcon'], programparams['samplefilter'],
                         programparams['varconfilter'], programparams['readidfilter'])
            
            # Run the LogInfo util.
            if utiltorun == "loginfo":
                li = LogInfo(self.vuh)
                li.main(programparams['vaselog'], programparams['logfilter'])

            # Subset a provided variant context file and write it to a new output file
            if utiltorun == "subsetvarcon":
                varcon_file = VariantContextFile(programparams["varcon"], programparams["samplefilter"],
                                                 programparams["varconfilter"], programparams["chromfilter"])
                varcon_file.write_variant_context_file(programparams["outfile"])
            
            # Run the VarconData util
            if utiltorun == "varcondata":
                varconvcfdata = VarconVcfData()
                varconvcfdata.main(programparams['donorfiles'], programparams['varcon'], programparams['samplefilter'],
                               programparams['varconfilter'], programparams['chromfilter'])
        else:
            self.vaseutillogger.warning("Not all parameters were set.")
            notsetparams = self.upc.get_not_set_parameters(utiltorun, programparams)
            self.vaseutillogger.warning("Parameter(s) " + ", ".join(notsetparams) + " are invalid")
            self.vaseutillogger.warning("Skipping selected util")

    # Start the VaSeUtil logger to be used for logging purposes
    def start_logger(self):
        vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        return vaseutillogger


# Run VaSeUtils
vsu = VaSeUtils()
vsu.main()
