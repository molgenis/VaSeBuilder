#!/usr/bin/env python
import logging
import argparse
import uuid
import sys

# Import the classes 
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext
from GvcfFile import GvcfFile
from GvcfVariant import GvcfVariant
from ParamChecker import ParamChecker
from VaSeEvaluator import VaSeEvaluator

class VaSeEval:
	# Constructor that only calls the main method.
	def __init__(self):
		self.main()
	
	
	# Performs the actual work of the program.
	def main(self):
		argList = self.getVaSeEvalParameters()
		pmc = ParamChecker()
		self.vaseEvalLogger = self.startLogger(pmc, argList['log'])
		
		if(pmc.check_parameters(argList)):
			self.vaseEvalLogger.info("Start NGSD_DNA pipeline result evaluation for datasets created by VaSeBuilder")
			vse = VaSeEvaluator(uuid.uuid4().hex)
			#vse.performVariantCallingEvaluation(argList)
			#self.vaseEvalLogger.info("VaSe evaluation completed")
			print("aap")
		else:
			print("jan")
			self.vaseEvalLogger.critical("Not all required parameters are ok. Check log for more info")
	
	
	# Returns the name and value of the command line arguments.
	def getVaSeEvalParameters(self):
		vaseEvalArgPars = argparse.ArgumentParser()
		vaseEvalArgPars.add_argument("--varcon", required=True, help="Location to write variants and their contexts to.", metavar="VARCON")
		vaseEvalArgPars.add_argument("--varbread", required=True, help="Location to write the variants and associated BAM reads to.", metavar="VARBREAD")
		vaseEvalArgPars.add_argument("--acceptorbread", required=True, help="Location to write the variants and associated template BAM reads to.", metavar="TEMPLATEBREAD")
		#vaseEvalArgPars.add_argument("--gvcf", required=True, help="Location of the gVCF file produced by the pipeline.", metavar="GVCF")
		#vaseEvalArgPars.add_argument("--donorvcfin", required=True, help="File containing a list of used donor VCF files.", metavar="VCFIN")
		#vaseEvalArgPars.add_argument("--donorbamin", required=True, help="File containing a list of used donor BAM files.", metavar="BAMIN")
		vaseEvalArgPars.add_argument("--out", required=True, help="Folder to write all output files to.", metavar="OUT")
		vaseEvalArgPars.add_argument("--log", help="Location to write log files to (will write to working directory if not used).", metavar="LOGFILE")
		return vars(vaseEvalArgPars.parse_args())
	
	
	# Check that all the parameters are ok.
	def parametersAreOk(self, argList):
		pmc = ParamChecker()
		return pmc.check_parameters(argList)
	
	
	# Method that creates the logger thagt will write the log to stdout and a log file.
	def startLogger(self, paramCheck, logloc):
		vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
		vaseEvalLogger.setLevel(logging.INFO)
		vaseLogFormat = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
		
		# Add the log stream to stdout
		vaseCliHandler = logging.StreamHandler(sys.stdout)
		vaseCliHandler.setLevel(logging.INFO)
		vaseCliHandler.setFormatter(vaseLogFormat)
		vaseEvalLogger.addHandler(vaseCliHandler)
		
		# Create the log stream to log file.
		logloc = paramCheck.check_log(logloc)
		if(logloc == ""):
			logloc = "VaSeEval.log"
		vaseFileHandler = logging.FileHandler(logloc)
		vaseFileHandler.setLevel(logging.INFO)
		vaseFileHandler.setFormatter(vaseLogFormat)
		vaseEvalLogger.addHandler(vaseFileHandler)
		return vaseEvalLogger

# Run the program.
VaSeEval()
