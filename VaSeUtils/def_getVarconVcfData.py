import pysam
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

# Reads the donor list file and saves it into a 
def readUsedDonorListFile(dListFileLoc):
	donorList = {}
	try:
		with open(dListFileLoc, 'r') as dlistFile:
			next(dlistFile)
			for fileLine in dlistFile:
				fileLine = fileLine.strip()
				fileLineData = fileLine.split()
				donorList[fileLineData[0]] = fileLineData[1:].split(' ; ')
				
	except IOError as ioe:
		print("Could not read donor list file.")
		exit()
	return donorList


# 
def getVarconVcfData(dListFileLoc, varconFileLoc, sampleFilter=None, chromFilter=None, posFilter=None):
	donorList = readUsedDonorListFile(dListFileLoc)
	varconFile = VariantContextFile(varconFileLoc, sampleFilter, chromFilter, posFilter)
	variantContexts = varconFile.get_variant_contexts()
	
	for dvcfFileSample in donorList:
		if(dvcfFileSample in sampleFilter):
			processVcfFile(dvcfFileSample, donorList[dvcfFileSample])


# Processes a single VCF file.
def processVcfFile(sampleId, varconData, vcfFileLoc):
	try:
		vcfFile = pysam.VariantFile(vcfFileLoc)
		for varcon in varconData:
			vcfFile.fetch(varcon.get_variant_context_chrom(), )
		vcfFile.close()
	except IOError as ioe:
		print("Could not read VCF file for sample " +str(vcfFileLoc))