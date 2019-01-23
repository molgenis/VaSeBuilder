#!/usr/bin/python

#Import necessary modules
import logging
import sys
import os
from datetime import datetime
import uuid
import argparse
from Bio import SeqIO
import gzip

#Import VaSe classes
from ParamChecker import ParamChecker
from VcfBamScanner import VcfBamScanner
from VaSeBuilder import VaSeBuilder


#Method that creates the logger thagt will write the log to stdout and a log file.
def startLogger(paramCheck, logloc):
	vaseLogger = logging.getLogger("VaSe_Logger")
	vaseLogger.setLevel(logging.INFO)
	vaseLogFormat = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
	
	#Add the log stream to stdout
	vaseCliHandler = logging.StreamHandler(sys.stdout)
	vaseCliHandler.setLevel(logging.INFO)
	vaseCliHandler.setFormatter(vaseLogFormat)
	vaseLogger.addHandler(vaseCliHandler)
	
	#Create the log stream to log file.
	logloc = paramCheck.checkLog(logloc)
	if(logloc == ""):
		logloc = "VaSeBuilder.log"
	vaseFileHandler = logging.FileHandler(logloc)
	vaseFileHandler.setLevel(logging.INFO)
	vaseFileHandler.setFormatter(vaseLogFormat)
	vaseLogger.addHandler(vaseFileHandler)


#Define the parameters and obtain their values.
vaseArgPars = argparse.ArgumentParser()
vaseArgPars.add_argument("--vcfin", nargs="*", required=True, help="Folder(s) containing VCF files to use.", metavar="VCFIN")
vaseArgPars.add_argument("--bamin", nargs="*", required=True, help="Folder(s) containing BAM files to use.", metavar="BAMIN")
vaseArgPars.add_argument("--bam", required=True, help="Location of the BAM file to modify and produce new FastQ.", metavar="BAM")
vaseArgPars.add_argument("--fastqout", required=True, help="Location to write the FastQ output file to.", metavar="FOUT")
vaseArgPars.add_argument("--varcon", help="Location to write variants and their contexts to.", metavar="VARCON")
vaseArgPars.add_argument("--varbread", help="Location to write the variants and associated BAM reads to.", metavar="VARBREAD")
vaseArgPars.add_argument("--log", help="Location to write log files to (will write to working directory if not used).", metavar="LOGFILE")
vaseArgs = vaseArgPars.parse_args()


#Make the required objects.
pmc = ParamChecker()
vaseLogger = startLogger(pmc, vars(vaseArgs)['log'])
vbscan = VcfBamScanner()
vaseB = VaSeBuilder(uuid.uuid4().hex)


#Proceed if the parameters are ok.
if(pmc.checkParameters(vars(vaseArgs))):
	#Start scanning the VCF and BAM folders
	vcfFiles = vbscan.scanVcfFolders(pmc.getValidVcfFolders())
	bamFiles = vbscan.scanBamFolders(pmc.getValidBamFolders())
	
	#Start the procedure to build the validation set.
	vaseB.buildValidationSet(vcfFiles, bamFiles, pmc.getNistBam(), pmc.getFastqOut())
	vaseLogger.info("VaSeBuilder run completed succesfully.")
else:
	vaseLogger.critical("Not all parameters are correct. Please check log for more info.")
	exit()
