#!/usr/bin/env python

# Import necessary modules
import io
import logging
import sys
import os
from datetime import datetime
import uuid
import argparse
#import gzip
#import pysam

# Import VaSe classes
from ParamChecker import ParamChecker
from VcfBamScanner import VcfBamScanner
from VaSeBuilder import VaSeBuilder


class VaSe:
	# Performs the check that VaSe is run with Python 3.x
	def __init__(self):
		if(sys.version_info[0] < 3):
			raise Exception("Please run this program in Python 3")
			exit()
	
	
	# Runs the VaSeBuilder program.
	def main(self):
		vaseArgList = self.getVaSeParameters()
		pmc = ParamChecker()
		self.vaseLogger = self.startLogger(pmc, vaseArgList['log'])
		
		if(pmc.checkParameters(vaseArgList)):
			vbscan = VcfBamScanner()
			vaseB = VaSeBuilder(uuid.uuid4().hex)
			
			# Start scanning the VCF and BAM folders
			vcfFileMap = vbscan.scanVcfFolders(pmc.getValidVcfFolders())
			bamFileMap = vbscan.scanBamFolders(pmc.getValidBamFolders())
			vcfBamFileLinker = vbscan.getVcfToBamMap()
			
			if(len(vcfBamFileLinker) > 0):
				# Start the procedure to build the validation set.
				vaseB.buildValidationSet(vcfBamFileLinker, vcfFileMap, bamFileMap, pmc.getAcceptorBam(), pmc.getFirstFastqInLocation(), pmc.getSecondFastqInLocation(), pmc.getOutDirLocation(), pmc.getFastqOutLocation(), pmc.getVariantContextOutLocation())
				self.vaseLogger.info("VaSeBuilder run completed succesfully.")
			else:
				self.vaseLogger.critical("No valid samples available to create new validation set")
		else:
			self.vaseLogger.critical("Not all parameters are correct. Please check log for more info.")
			exit()
	
	
	# Method that creates the logger thagt will write the log to stdout and a log file.
	def startLogger(self, paramCheck, logloc, debug=False):
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
	
	
	# Returns the vase Parameters
	def getVaSeParameters(self):
		# Set the VaSe parameters for the program
		vaseArgPars = argparse.ArgumentParser()
		vaseArgPars.add_argument("-v", "--donorvcf", dest="donorvcf", nargs="*", required=True, help="Folder(s) containing VCF files.", metavar="DONORVCF")
		vaseArgPars.add_argument("-b", "--donorbam", dest="donorbam", nargs="*", required=True, help="Folder(s) containing BAM files.", metavar="DONORBAM")
		vaseArgPars.add_argument("-a", "--acceptorbam", dest="acceptorbam", required=True, help="Location of the BAM file to modify and produce new FastQ.", metavar="ACCEPTORBABM")
		vaseArgPars.add_argument("-1", "--templatefq1", dest="templatefq1", nargs="*", required=True, help="Location and name of the first fastq in file.", metavar="TEMPLATEFQ1")
		vaseArgPars.add_argument("-2", "--templatefq2", dest="templatefq2", nargs="*", required=True, help="Location and name of the second fastq in file.", metavar="TEMPLATEFQ2")
		vaseArgPars.add_argument("-o", "--out", dest="out", required=True, help="Location to write the output files to", metavar="OUT")
		vaseArgPars.add_argument("-of", "--fastqout", dest="fastqout", help="Name for the two FastQ files to be produced.", metavar="FASTQOUT")
		vaseArgPars.add_argument("-ov", "--varcon", dest="varcon", help="File name to write variants and their contexts to.", metavar="VARCON")
		vaseArgPars.add_argument("-l", "--log", dest="log", help="Location to write log files to (will write to working directory if not used).", metavar="LOGFILE")
		vaseArgPars.add_argument("-!", "--debug", dest='debug', default=False, action="store_true", help="Run the program in debug mode", metavar="DEBUG")
		vaseArgs = vars(vaseArgPars.parse_args())
		return vaseArgs

# Run the program	
vaseb = VaSe()
vaseb.main()