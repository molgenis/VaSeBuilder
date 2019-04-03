import argparse

from AcceptorCheck import AcceptorCheck
from AcceptorCheck import DonorCheck

class VaSeUtils:
	def main():
		vaseuArgs = self.getVaSeUtilsParameters()
		# Should check the parameters here!
		self.vaseUtilLogger.info("Running selected VaSeUtil program(s)")
		for utilToRun in vaseuArgs['util']:
			self.runSelectedUtil(utilToRun, vaseuArgs)
		self.vaseUtilLogger.info("Ran selected VaSeUtil program(s)")
	
	
	# Returns the set parameter values.
	def getVaSeUtilsParameters(self):
		vaseuArgPars = argparse.ArgumentParser(description="Run a specific VaSe Util program")
		vaseuArgPars.add_argument("-u", "--util", nargs="*", choices=['acceptorcheck', 'acceptorreadinfo', 'checkfastq', 'donorcheck', 'acceptorreadinfo', 'varcondata'], required=True, help="The utility to run.", metavar="UTIL")
		vaseuArgpars.add_argument("-df", "--donorfiles", help="File containing the used donor VCF/BAM files", metavar="DONORFILES")
		vaseuArgpars.add_argument("-dr", "--donorreads", help="File containing the donor reads per variant context", metavar="DONORREADS")
		vaseuArgpars.add_argument("-ar", "--acceptorreads", help="File containing the acceptor reads per variant context", metavar="ACCEPTORREADS")
		vaseuArgpars.add_argument("-vf1", "--vasefastq1", help="The VaSe produced R1 FastQ file", metavar="VASEFASTQ1")
		vaseuArgpars.add_argument("-vf2", "--vasefastq2", help="The VaSe produced R2 FastQ file", metavar="VASEFASTQ2")
		vaseuArgpars.add_argument("-tf1", "--templatefastq1", help="Template R1 FastQ file used to produce the VaSe R1 FastQ file", metavar="TEMPLATEFASTQ1")
		vaseuArgpars.add_argument("-tf2", "--templatefastq2", help="Template R2 FastQ file used to produce the VaSe R2 FastQ file", metavar="TEMPLATEFASTQ2")
		vaseuArgpars.add_argument("-ab", "--acceptorbam", help="BAM file used as acceptor", metavar="ACCEPTORBAM")
		vaseuArgpars.add_argument("-sf", "--samplefilter", help="List of sample identifiers to include. Will use all samples if not set", metavar="SAMPLEFILTER")
		vaseuArgpars.add_argument("-cf", "--chromfilter", help="List of chromosomes to use. Will use all chromosomes if not set", metavar="CHROMFILTER")
		vaseuArgpars.add_argument("-pf", "--posfilter", help="List of start|end positions to use. Will use all positions if not set", metavar="POSFILTER")
		vaseuArgPars.add_argument("-vf", "--varconfilter", help="List of variant context to use. Will use all variant contexts if not set", metavar="VARCONFILTER")
		return vars(vaseuArgPars.parse_args())
	
	
	# Runs one or more selected programs
	def runSelectedUtil(utilToRun, programParams):
		if(utilToRun=='acceptorcheck'):
			acheck = AcceptorCheck()
			bamReadList = self.readBamReadsList_noFilter(programParams['acceptorreads'])
			acheck.main(bamReadList, programParams['vasefastq1'], programParams['vasefastq2'])
		
		if(acceptorreadinfo=='acceptorreadinfo'):
			adri = AdReadInfo()
			adri.getReadInfo()
		
		if(utilToRun=='checkfastq'):
			acceptorReadList = self.readBamReadsList_noFilter(programParams['acceptorreads'])
			donorReadList = self.readBamReadsList_noFilter(programParams['donorreads'])
			checkf = CheckVaSeFastq()
			checkf.main()
		
		if(utilToRun=='donorcheck'):
			dcheck = DonorCheck()
			bamReadList = self.readBamReadsList_noFilter(programParams['donorreads'])
			dcheck.main(bamReadList, programParams['vasefastq1'], programParams['vasefastq2'])
		
		if(utilToRun=='donorreadinfo'):
			adri = AdReadInfo()
			adri.getReadInfo()
		
		if(utilToRun=='varcondata'):
			print("aap")
	
	
	# Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck	').
	def readBamReadsList_noFilter(self, acceptorReadFile):
		acceptorReads = []
		with open(acceptorReadFile, 'r') as arFile:
			next(arFile)	# Skip the header line
			for fileLine in arFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split("\t")
				acceptorReads.extend(fileLineData[1:])
		return acceptorReads


# Run VaSeUtils
vsu = VaSeUitls()
vsu.main()