#!/usr/bin/env python
import logging
import argparse

# Import the classes 
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext
from GvcfFile import GvcfFile
from GvcfVariant import GvcfVariant
from ParamChecker import ParamChecker

class VaSeEval:
	# Constructor that only calls the main method.
	def __init__(self):
		self.__main__()
	
	
	# Performs the actual work of the program.
	def __main__(self):
		argList = self.getVaSeEvalParameters()
		if(self.parametersAreOk(argList)):
			self.performValidation()	# Do the validation work
		else:
			self.vaseEvalLogger.critical("Not all required parameters are OK.")
			exit()
	
	
	# Returns the name and value of the command line arguments.
	def getVaSeEvalParameters(self):
		vaseEvalArgPars = argparse.ArgumentParser()
		vaseEvalArgPars.add_argument("--varcon", required=True, help="Location to write variants and their contexts to.", metavar="VARCON")
		vaseEvalArgPars.add_argument("--varbread", required=True, help="Location to write the variants and associated BAM reads to.", metavar="VARBREAD")
		vaseEvalArgPars.add_argument("--acceptorbread", required=True, help="Location to write the variants and associated template BAM reads to.", metavar="TEMPLATEBREAD")
		vaseEvalArgPars.add_argument("--gvcf", required=True, help="Location of the gVCF file produced by the pipeline.", metavar="GVCF")
		vaseEvalArgPars.add_argument("--donorvcfin", required=True, help="File containing a list of used donor VCF files.", metavar="VCFIN")
		vaseEvalArgPars.add_argument("--donorbamin", required=True, help="File containing a list of used donor BAM files.", metavar="BAMIN")
		vaseEvalArgPars.add_argument("--out", required=True, help="Folder to write all output files to.", metavar="OUT")
		vaseEvalArgPars.add_argument("--log", help="Location to write log files to (will write to working directory if not used).", metavar="LOGFILE")
		return vars(vaseEvalArgPars.parse_args())
	
	
	# Check that all the parameters are ok.
	def parametersAreOk(self, argList):
		pmc = ParamChecker()
		return pmc.checkParameters(argList)
	
	
	# Perform the validtion work
	def performValidation(self, argList):
		# 1.Read the varcon, varbread and templatebread files
		varConFile = VariantContextFile(argList['varcon'])
		
		# 2.Read the gvcf file containing the results of
		resultsGvcf = GvcfFile(argList['gvcf'])
		
		# 3.Iterate over the donor vcf files and check all donor variants
		donorVcfFiles = self.readDonorVcfIn(argList['donorvcfin'])
		variantCallingResults = self.checkDonorVariants(donorVcfFiles, resultsGvcf, varConFile)
		
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'cc', argList['out']+"/varcal_cc.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'ci', argList['out']+"/varcal_ci.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'nc', argList['out']+"/varcal_nc.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'ic', argList['out']+"/varcal_ic.txt")
		
		# 4.Read the varbread and templatebread files and add the read identifiers to the variant contexts.
		#self.addReadIdsToContexts(argList['varbread'], varConFile, 'donor')
		#self.addReadIdsToContexts(argList['acceptorbread'], varConFile, 'acceptor')
		
		# 5.Iterate over the donor bam files and check which reads are also in the bam produced by the pipeline
		#donorBamFiles = self.readDonorBamIn(argList['donorbamin'])
	
	
	# Checks whether the donor variants added by VaSe are found and are found correctly (0/0, 0/1, etc)
	def checkDonorVariants(self, dvcfListFile, gvcf, varconfile):
		dVcfCalling = {}	# Will save the pipeline calling results based on cc, ci, nc, ic
		for sampleid, dvcfFile in donorVcfFiles:
			dVcfCalling[sampleid] = dvcfFile.checkVariantsInOther(sampleid, gvcf, varconfile)
		return dVcfCalling
	
	
	# Reads the file containing the locations of donor VCF files used by the VaSe to construct the validation set.
	def readDonorVcfIn(self, donorVcfLoc):
		dvcfFiles = {}
		try:
			with open(donorVcfLoc, 'r') as dVcfFile:
				for fileLine in dVcfFile.readlines():
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					try:
						dvcfFiles[fileLineData[0]] = GvcfFile(fileLineData[1])
					except IOError as ioe:
						self.vaseEvalLogger.warning("Could not read VCF file " +str(fileLineData[1]))
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not read the file with donor vcf locations. Maybe it has been moved/deleted?")
			exit()
		if(len(dvcfFiles)==0):
			self.vaseEvalLogger.critical("No valid donor VCF files could be read")
			exit()
		return dvcfFiles
	
	
	# Reads the file containing the locations of the donor BAM files used by VaSe to construct the validation set.
	def readDonorBamIn(self donorBamLoc):
		dbamFiles = {}
		try:
			with open() as dBamFile:
				for fileLine in dBamFile.readlines():
					fileLine = fileLine.strip()
					fileLineData = fileLine.split("\t")
					dbamFiles[fileLineData[0]] = fileLineData[1]
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not read the file with donor bam locations. Maybe it has been moved/deleted?")
			exit()
		
	
	# Adds donor/acceptor read identifiers to the variant contexts
	def addReadIdsToContexts(fileLoc, vcFile, donacc):
		try:
			with open(fileLoc, 'r') as varconReadFile:
				for fileLine in varconReadFile:
					fileLine = fileLine.strip()
					
					if(not fileLine.startswith("Variant")):
						fileLineData = fileLine.split("\t")
						if(donacc=='donor'):
							vcFile.addDonorReadsToVariantContext(fileLineData[0], fileLineData[1:])
						else:
							vcFile.addAcceptorReadIdentifiers(fileLineData[0], fileLineData[1:])
		except IOError as ioe:
			self.vaseEvalLogger.critical("Could not read bread file " +str(fileLoc))
			exit()
	
	
	# Writes variant calling results of cc, ci, nc and ic to a specified output file.
	def writeVariantCallingResults(self, varcalResults, donorVcfFiles, caltype, varcalOutLoc):
		try:
			with open(varcalOutLoc, 'w') as varcalFile:
				varcalFile.write("VariantId\tDonor_Ref\tCall_Ref\tDonor_Alt\tCall_Alt\tDonor_GT\tCall_GT\n")	# Write an understandable header line.
				for sampleId in varcalResults:
					variantList = varcalResults[sampleId][caltype]
					
					# Iterate over the variants and their calling result.
					for variantId in variantList:
						donorRef = donorVcfFiles[sampleId].getVariantRef(variantId)
						donorAlt = donorVcfFiles[sampleId].getVariantAlt(variantId)
						donorGt = donorVcfFiles[sampleId].getSampleInfo('GT')
						
						# Check if the caltype is not 'nc' or 'Not Called'.
						if(caltype != 'nc'):
							callRef = donorVcfFiles[sampleId].getCalledRef()
							callAlt = donorVcfFiles[sampleId].getCalledAlt()
							callGt = donorVcfFiles[sampleId].getCalledGt()
						else:
							callRef, callAlt, callGt = 'NA'
						varcalFile.write(str(variantId)+ "\t" +str(donorRef)+ "\t" +str(callRef)+ "\t" +str(donorAlt)+ "\t" +str(callAlt)+ "\t" +str(donorGt)+ "\t" +str(callGt)+ "\n")
		except IOerror as ioe:
			self.vaseEvalLogger.critical("Could not write calling results to")
			exit()


# Run the program.
VaSeEval()
