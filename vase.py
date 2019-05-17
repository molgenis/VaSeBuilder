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
        if (sys.version_info[0] < 3 and sys.version_info[1] < 6):
            raise Exception("Please run this program in at least Python 3.6")
            exit()
        if (int(pysam.version.__version__.split('.')[0]) < 1
           and int(pysam.version.__version__.split('.')[1]) < 15):
            raise Exception("Please run this program with at least Pysam 0.15")
            exit()

    # Runs the VaSeBuilder program.
    def main(self):
        vaseArgList = self.getVaSeParameters()
        pmc = ParamChecker()
        self.vaseLogger = self.startLogger(pmc,
                                           vaseArgList['log'],
                                           vaseArgList['debug'])

        if (pmc.checkParameters(vaseArgList)):
            vbscan = VcfBamScanner()
            vaseB = VaSeBuilder(uuid.uuid4().hex)

            # Start scanning the VCF and BAM folders.
            vcfFileMap = vbscan.scanVcfFolders(pmc.getValidVcfFolders())
            bamFileMap = vbscan.scanBamFolders(pmc.getValidBamFolders())
            vcfBamFileLinker = vbscan.getVcfToBamMap()

            if (len(vcfBamFileLinker) > 0):
                # Start the procedure to build the validation set.
                vaseB.buildValidationSet(vcfBamFileLinker,
                                         vcfFileMap, bamFileMap,
                                         pmc.getAcceptorBam(),
                                         pmc.getFirstFastqInLocation(),
                                         pmc.getSecondFastqInLocation(),
                                         pmc.getOutDirLocation(),
                                         pmc.getFastqOutLocation(),
                                         pmc.getVariantContextOutLocation())
                self.vaseLogger.info("VaSeBuilder run completed succesfully.")
            else:
                self.vaseLogger.critical("No valid samples available to "
                                         "create new validation set")
        else:
            self.vaseLogger.critical("Not all parameters are correct. Please "
                                     "check log for more info.")
            exit()

    # Method that creates the logger thagt will write the log to stdout
    # and a log file.
    def startLogger(self, paramCheck, logloc, debugMode=False):
        vaseLogger = logging.getLogger("VaSe_Logger")
        if (debugMode):
            vaseLogger.setLevel(logging.DEBUG)
        else:
            vaseLogger.setLevel(logging.INFO)
        vaseLogFormat = logging.Formatter(
                "%(asctime)s %(name)s %(levelname)s %(message)s"
                )

        # Add the log stream to stdout.
        vaseCliHandler = logging.StreamHandler(sys.stdout)
        if (debugMode):
            vaseCliHandler.setLevel(logging.DEBUG)
        else:
            vaseCliHandler.setLevel(logging.INFO)
        vaseCliHandler.setFormatter(vaseLogFormat)
        vaseLogger.addHandler(vaseCliHandler)

        # Create the log stream to log file.
        logloc = paramCheck.checkLog(logloc)
        if (logloc == ""):
            logloc = "VaSeBuilder.log"
        vaseFileHandler = logging.FileHandler(logloc)

        if (debugMode):
            vaseFileHandler.setLevel(logging.DEBUG)
        else:
            vaseFileHandler.setLevel(logging.INFO)
        vaseFileHandler.setFormatter(vaseLogFormat)
        vaseLogger.addHandler(vaseFileHandler)
        return vaseLogger

    # Returns the vase Parameters.
    def getVaSeParameters(self):
        # Set the VaSe parameters for the program.
        vaseArgPars = argparse.ArgumentParser()
        vaseArgPars.add_argument("-v", "--donorvcf", dest="donorvcf", nargs="+", required=True, help="Folder(s) containing VCF files.")
        vaseArgPars.add_argument("-b", "--donorbam", dest="donorbam", nargs="+", required=True, help="Folder(s) containing BAM files.")
        vaseArgPars.add_argument("-a", "--acceptorbam", dest="acceptorbam", required=True, help="Location of the BAM file to modify and produce new FastQ.")
        vaseArgPars.add_argument("-1", "--templatefq1", dest="templatefq1", nargs="+", required=True, help="Location and name of the first fastq in file.")
        vaseArgPars.add_argument("-2", "--templatefq2", dest="templatefq2", nargs="+", required=True, help="Location and name of the second fastq in file.")
        vaseArgPars.add_argument("-o", "--out", dest="out", required=True, help="Location to write the output files to")
        vaseArgPars.add_argument("-of", "--fastqout", dest="fastqout", help="Name for the two FastQ files to be produced.")
        vaseArgPars.add_argument("-ov", "--varcon", dest="varcon", help="File name to write variants and their contexts to.")
        vaseArgPars.add_argument("-l", "--log", dest="log", help="Location to write log files to (will write to working directory if not used).")
        vaseArgPars.add_argument("-!", "--debug", dest='debug', default=False, action="store_true", help="Run the program in debug mode")
        vaseArgs = vars(vaseArgPars.parse_args())
        return vaseArgs


# Run the program.
vaseb = VaSe()
vaseb.main()
