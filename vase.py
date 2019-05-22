#!/usr/bin/env python

# Import necessary modules.
import io
import logging
import sys
import os
from datetime import datetime
import uuid
import argparse
import gzip
import pysam

# Import VaSe classes.
from ParamChecker import ParamChecker
from VcfBamScanner import VcfBamScanner
from VaSeBuilder import VaSeBuilder


class VaSe:
    # Performs the check that VaSe is run with Python 3.x
    def __init__(self):
        assert (sys.version_info[0] >= 3 and sys.version_info[1] >= 6), "Please run this program in Python 3.6 or higher"
        assert (int(pysam.version.__version__.split('.')[0]) >= 0 and int(pysam.version.__version__.split('.')[1]) >= 15), "Please run this program with Pysam 0.15 or higher"

    # Runs the VaSeBuilder program.
    def main(self):
        vase_arg_list = self.get_vase_parameters()
        pmc = ParamChecker()
        self.vaselogger = self.start_logger(pmc,
                                            vase_arg_list['log'],
                                            vase_arg_list['debug'])

        if pmc.check_parameters(vase_arg_list):
            vbscan = VcfBamScanner()
            vase_b = VaSeBuilder(uuid.uuid4().hex)

            # Start scanning the VCF and BAM folders.
            vcf_file_map = vbscan.scan_vcf_folders(pmc.get_valid_vcf_folders())
            bam_file_map = vbscan.scan_bam_folders(pmc.get_valid_bam_folders())
            vcf_bam_file_linker = vbscan.get_vcf_to_bam_map()

            if len(vcf_bam_file_linker) > 0:
                # Start the procedure to build the validation set.
                vase_b.build_validation_set(vcf_bam_file_linker,
                                            vcf_file_map, bam_file_map,
                                            pmc.get_acceptor_bam(),
                                            pmc.get_first_fastq_in_location(),
                                            pmc.get_second_fastq_in_location(),
                                            pmc.get_out_dir_location(),
                                            pmc.get_fastq_out_location(),
                                            pmc.get_variant_context_out_location())
                self.vaselogger.info("VaSeBuilder run completed succesfully.")
            else:
                self.vaselogger.critical("No valid samples available to "
                                         "create new validation set")
        else:
            self.vaselogger.critical("Not all parameters are correct. Please "
                                     "check log for more info.")
            exit()

    # Method that creates the logger thagt will write the log to stdout
    # and a log file.
    def start_logger(self, paramcheck, logloc, debug_mode=False):
        vaselogger = logging.getLogger("VaSe_Logger")
        if debug_mode:
            vaselogger.setLevel(logging.DEBUG)
        else:
            vaselogger.setLevel(logging.INFO)
        vaselog_format = logging.Formatter(
                "%(asctime)s %(name)s %(levelname)s %(message)s"
                )

        # Add the log stream to stdout.
        vase_cli_handler = logging.StreamHandler(sys.stdout)
        if debug_mode:
            vase_cli_handler.setLevel(logging.DEBUG)
        else:
            vase_cli_handler.setLevel(logging.INFO)
        vase_cli_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_cli_handler)

        # Create the log stream to log file.
        logloc = paramcheck.check_log(logloc)
        if logloc == "":
            logloc = "VaSeBuilder.log"
        vase_file_handler = logging.FileHandler(logloc)

        if debug_mode:
            vase_file_handler.setLevel(logging.DEBUG)
        else:
            vase_file_handler.setLevel(logging.INFO)
        vase_file_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_file_handler)
        return vaselogger

    # Returns the vase Parameters.
    def get_vase_parameters(self):
        # Set the VaSe parameters for the program.
        vase_argpars = argparse.ArgumentParser()
        vase_argpars.add_argument("-v", "--donorvcf", dest="donorvcf", nargs="+", required=True, help="Folder(s) containing VCF files.")
        vase_argpars.add_argument("-b", "--donorbam", dest="donorbam", nargs="+", required=True, help="Folder(s) containing BAM files.")
        vase_argpars.add_argument("-a", "--acceptorbam", dest="acceptorbam", required=True, help="Location of the BAM file to modify and produce new FastQ.")
        vase_argpars.add_argument("-1", "--templatefq1", dest="templatefq1", nargs="+", required=True, help="Location and name of the first fastq in file.")
        vase_argpars.add_argument("-2", "--templatefq2", dest="templatefq2", nargs="+", required=True, help="Location and name of the second fastq in file.")
        vase_argpars.add_argument("-o", "--out", dest="out", required=True, help="Location to write the output files to")
        vase_argpars.add_argument("-of", "--fastqout", dest="fastqout", help="Name for the two FastQ files to be produced.")
        vase_argpars.add_argument("-ov", "--varcon", dest="varcon", help="File name to write variants and their contexts to.")
        vase_argpars.add_argument("-l", "--log", dest="log", help="Location to write log files to (will write to working directory if not used).")
        vase_argpars.add_argument("-!", "--debug", dest='debug', default=False, action="store_true", help="Run the program in debug mode")
        vase_args = vars(vase_argpars.parse_args())
        return vase_args


# Run the program.
vaseb = VaSe()
vaseb.main()
