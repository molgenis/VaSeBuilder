import argparse

from AcceptorCheck import AcceptorCheck
from checkVaseFastQ import CheckVaSeFastQ
from donorcheck import DonorCheck
from DonorReadInfo import DonorReadInfo
from UtilParamCheck import UtilParamCheck
from VaSeUtilHelper import VaSeUtilHelper

class VaSeUtils:
	def __init__(self):
		self.upc = UtilParamCheck()
		self.vuh = VaSeUtilHelper()
		self.validUtils = ['acceptorcheck', 'acceptorreadinfo', 'checkfastq', 'compareacceptor', 'comparedonor', 'comparefastq', 'comparevarcon', 'donorcheck', 'donorreadinfo', 'loginfo', 'unmappedinfo', 'varcondata', 'vasecompare']
	
	# Runs all specified VaSe utils
	def main():
		vaseuArgs = self.getVaSeUtilsParameters()
		# Should check the parameters here!
		self.vaseUtilLogger.info("Running selected VaSeUtil program(s)")
		
		# Run one or more sleected utils
		for utilToRun in vaseuArgs['util']:
			self.runSelectedUtil(utilToRun, vaseuArgs)
		self.vaseUtilLogger.info("Ran selected VaSeUtil program(s)")
	
	
	# Returns the set parameter values.
	def getVaSeUtilsParameters(self):
		vaseuArgPars = argparse.ArgumentParser(description="Run a specific VaSe Util program")
		vaseuArgPars.add_argument("-u", "--util", nargs="*", dest='util', choices=self.validUtils, required=True, help="The utility to run.", metavar="UTIL")
		vaseuArgPars.add_argument("-l", "--log", dest='log', help="Location to write log file to", metavar="LOGFILE")
		vaseuArgPars.add_argument("-df", "--donorfiles", dest='donorfiles', help="File containing the list of used donor VCF/BAM files", metavar="DONORFILES")
		vaseuArgPars.add_argument("-dr", "--donorreads", dest='donorreads', help="File containing the donor reads per variant context", metavar="DONORREADS")
		vaseuArgPars.add_argument("-ar", "--acceptorreads", dest='acceptorreads', help="File containing the acceptor reads per variant context", metavar="ACCEPTORREADS")
		vaseuArgPars.add_argument("-vf1", "--vasefastq1", dest='vasefastq1', help="The VaSe produced R1 FastQ file", metavar="VASEFASTQ1")
		vaseuArgPars.add_argument("-vf2", "--vasefastq2", dest='vasefastq2', help="The VaSe produced R2 FastQ file", metavar="VASEFASTQ2")
		vaseuArgPars.add_argument("-vl", "--vaselog", dest='vaselog', help="Location to the log file produced by VaSeBuilder", metavar="VASELOG")
		vaseuArgPars.add_argument("-tf1", "--templatefastq1", dest='templatefastq1', help="Template R1 FastQ file used to produce the VaSe R1 FastQ file", metavar="TEMPLATEFASTQ1")
		vaseuArgPars.add_argument("-tf2", "--templatefastq2", dest='templatefastq2', help="Template R2 FastQ file used to produce the VaSe R2 FastQ file", metavar="TEMPLATEFASTQ2")
		vaseuArgPars.add_argument("-ab", "--acceptorbam", dest='acceptorbam', help="BAM file used as acceptor", metavar="ACCEPTORBAM")
		vaseuArgPars.add_argument("-vc", "--varcon", dest='varcon', help="VaSe produced variant context file", metavar="VARCON")
		vaseuArgPars.add_argument("-vc2", "--varcon2", dest='varcon2', help="Other VaSe produced variant context file", metavar="VARCON2")
		vaseuArgPars.add_argument("-um", "--unmappedmates", dest='unmappedmates', help="VaSe produced file with reads that have unmapped mates", metavar="UNMAPPEDMATES")
		vaseuArgPars.add_argument("-um2", "--unmappedmates2", dest='unmappedmates2', help="Other VaSe produced file with reads that have unmapped mates", metavar="UNMAPPEDMATES2")
		vaseuArgPars.add_argument("-sf", "--samplefilter", dest='samplefilter', help="List of sample identifiers to include. Will use all samples if not set", metavar="SAMPLEFILTER")
		vaseuArgPars.add_argument("-cf", "--chromfilter", dest='chromfilter', help="List of chromosomes to use. Will use all chromosomes if not set", metavar="CHROMFILTER")
		vaseuArgPars.add_argument("-pf", "--posfilter", dest='posfilter', help="List of start-end position ranges to use. Will use all positions if not set", metavar="POSFILTER")
		vaseuArgPars.add_argument("-vf", "--varconfilter", dest='varconfilter', help="List of variant context to use. Will use all variant contexts if not set", metavar="VARCONFILTER")
		vaseuArgpars.add_argument("-lf", "--logfilter", dest='logfilter', help="Filter for which log fields to show (e.g. INFO, DEBUG, WARNING)", metavar="LOGFILTER")
		return vars(vaseuArgPars.parse_args())
	
	
	# Runs a selected util
	def runSelectedUtil(utilToRun, programParams):
		if(upc.requiredParamsSet(utilToRun, programParams)):
			if(utilToRun=='acceptorcheck'):
				acheck = AcceptorCheck()
				bamReadList = self.vuh.readBamReadsList_noFilter(programParams['acceptorreads'])
				acheck.main(bamReadList, programParams['vasefastq1'], programParams['vasefastq2'])
			
			if(acceptorreadinfo=='acceptorreadinfo'):
				ari = AcceptorReadInfo(self.vuh)
				ari.main(programParams['acceptorbam'], programParams['acceptorreads'], programParams['varcon'], programParams['samplefilter'], programParams['varconfilter'])
			
			if(utilToRun=='checkfastq'):
				acceptorReadList = self.vuh.readABamReadsList_noFilter(programParams['acceptorreads'])
				donorReadList = self.vuh.readDBamReadsList_noFilter(programParams['donorreads'])
				checkf = CheckVaSeFastq()
				checkf.main()
			
			if(utilToRun=='donorcheck'):
				dcheck = DonorCheck()
				bamReadList = self.vuh.readBamReadsList_noFilter(programParams['donorreads'])
				dcheck.main(bamReadList, programParams['vasefastq1'], programParams['vasefastq2'])
			
			if(utilToRun=='donorreadinfo'):
				dri = DonorReadInfo(self.vuh)
				dri.main()
			
			if(utilToRun=='loginfo'):
				print("aap")
			
			if(utilToRun=='varcondata'):
				print("aap")
		else:
			self.vaseUtilLogger.warning("Not all parameters were set.")
			notSetParams = self.upc.getNotSetParameters(utilToRun, programParams)
			self.vaseUtilLogger.warning("Parameter(s) " +", ".join(notSetParams)+ " are invalid")

# Run VaSeUtils
vsu = VaSeUtils()
vsu.main()