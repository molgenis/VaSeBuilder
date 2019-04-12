#!/usr/bin/env python
import logging
from datetime import datetime

class VaSeEvaluator:
	
	# Set the identifier and creation time and date of this vase evaluator
	def __init__(self, vseId):
		self.vaseEvalId = vseId
		self.creationDate = datetime.now().date()
		self.creationTime = datetime.now().time()
		self.vaseEvalLogger = logging.getLogger("VaSeEval_Logger")
	
	
	# Performs the validation of the variant calling done by the pipeline
	def performVariantCallingEvaluation(self, argList):
		varConFile = VariantContextFile(argList['varcon'])	# 1.Read the varcon, varbread and templatebread files
		resultsGvcf = GvcfFile(argList['gvcf'])	# 2.Read the gvcf file containing the variant calling results produced by the pipeline
		
		donorVcfFiles = self.readDonorVcfIn(argList['donorvcfin'])
		variantCallingResults = self.checkDonorVariants(donorVcfFiles, resultsGvcf, varConFile)	# 3.Iterate over the donor vcf files and check all donor variants
		
		# Write the variant calling evaluations to four output files.
		self.vaseEvalLogger.info("Start writing the calling results for all variants")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'cc', argList['out']+"/varcal_cc.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'ci', argList['out']+"/varcal_ci.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'nc', argList['out']+"/varcal_nc.txt")
		self.writeVariantCallingResults(variantCallingResults, donorVcfFiles, 'ic', argList['out']+"/varcal_ic.txt")
	
	
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
	def readDonorBamIn(self, donorBamLoc):
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
	
	
	#def performBamReadValidation(self):
