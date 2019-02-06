#!/usr/bin/env python

# Import necessary modules
import logging
import sys
import os
from datetime import datetime
import uuid
import argparse
from Bio import SeqIO
import gzip
import pysam

# Import VaSe classes
from ParamChecker import ParamChecker
from VcfBamScanner import VcfBamScanner
from VaSeBuilder import VaSeBuilder


# Method that creates the logger thagt will write the log to stdout and a log file.
def startLogger(paramCheck, logloc):
	vaseLogger = logging.getLogger("VaSe_Logger")
	vaseLogger.setLevel(logging.INFO)
	vaseLogFormat = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
	
	# Add the log stream to stdout
	vaseCliHandler = logging.StreamHandler(sys.stdout)
	vaseCliHandler.setLevel(logging.INFO)
	vaseCliHandler.setFormatter(vaseLogFormat)
	vaseLogger.addHandler(vaseCliHandler)
	
	# Create the log stream to log file.
	logloc = paramCheck.checkLog(logloc)
	if(logloc == ""):
		logloc = "VaSeBuilder.log"
	vaseFileHandler = logging.FileHandler(logloc)
	vaseFileHandler.setLevel(logging.INFO)
	vaseFileHandler.setFormatter(vaseLogFormat)
	vaseLogger.addHandler(vaseFileHandler)
	return vaseLogger



# Define the parameters and obtain their values. (Maybe combine '--varcon', '--varbread' and '--nistbread' into one parameter such as '--out')
vaseArgPars = argparse.ArgumentParser()
vaseArgPars.add_argument("--vcfin", nargs="*", required=True, help="Folder(s) containing VCF files.", metavar="VCFIN")
vaseArgPars.add_argument("--bamin", nargs="*", required=True, help="Folder(s) containing BAM files.", metavar="BAMIN")
vaseArgPars.add_argument("--valbam", required=True, help="Location of the BAM file to modify and produce new FastQ.", metavar="VALBAM")
vaseArgPars.add_argument("--valfastq1", help="Location and name of the first fastq in file.", metavar="VALFIN1")
vaseArgPars.add_argument("--valfastq2", help="Location and name of the second fastq in file.", metavar="VALFIN2")
vaseArgPars.add_argument("--fastqout", required=True, help="Folder to write the new FastQ output file to.", metavar="FOUT")
vaseArgPars.add_argument("--varcon", required=True, help="Location to write variants and their contexts to.", metavar="VARCON")
vaseArgPars.add_argument("--varbread", required=True, help="Location to write the variants and associated BAM reads to.", metavar="VARBREAD")
vaseArgPars.add_argument("--nistbread", required=True, help="")
vaseArgPars.add_argument("--log", help="Location to write log files to (will write to working directory if not used).", metavar="LOGFILE")
vaseArgs = vaseArgPars.parse_args()


# Make the required objects.
pmc = ParamChecker()
vaseLogger = startLogger(pmc, vars(vaseArgs)['log'])
vbscan = VcfBamScanner()
vaseB = VaSeBuilder(uuid.uuid4().hex)


# Proceed if the parameters are ok.
if(pmc.checkParameters(vars(vaseArgs))):
	# Start scanning the VCF and BAM folders
	vcfFiles = vbscan.scanVcfFolders(pmc.getValidVcfFolders())
	bamFiles = vbscan.scanBamFolders(pmc.getValidBamFolders())
	
	vcfBamFileLinker = vbscan.getVcfToBamMap()
	
	if(len(vcfBamFileLinker) > 0):
		# Start the procedure to build the validation set.
		vaseB.buildValidationSet(vcfFileMap, bamFileMap, pmc.getNistBam(), pmc.getFirstFastqInLocation(), pmc.getSecondFastqInLocation(), pmc.getFastqOutLocation(), pmc.getVariantContextOutLocation(), pmc.getVariantBamReadOutLocation(), pmc.getNistBamReadOutLocation())
		vaseLogger.info("VaSeBuilder run completed succesfully.")
	else:
		vaseLogger.critical("No valid samples available to create new validation set")
else:
	vaseLogger.critical("Not all parameters are correct. Please check log for more info.")
	exit()
