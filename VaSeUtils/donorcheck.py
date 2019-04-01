import sys
import gzip

# Reads the list of added donor reads
def readDonorReadList(donorReadFile):
	donorReadList = []
	with open(donorReadFile, 'r') as drFile:
		next(drFile)	# Skip the header line of the file
		for fileLine in drFile:
			fileLine = fileLine.strip()
			fileLineData = fileLine.split("\t")
			donorReadList.extend(fileLineData[2].split(" ; "))
	return donorReadList


# Checks whether the list of reads are in a specified fastq file based on read identifier.
def checkDonorReadsAdded(gzResultsFile, donorReadFile):
	donorReadList = readDonorReadList(donorReadFile)
	addedCount = 0
	
	# Checks whether the identified donor reads have been added.
	with gzip.open(gzResultsFile, 'rt') as gzFile:
		for fileLine in gzFile:
			fileLine = fileLine.strip()
			if(fileLine.startswith('@')):
				if(fileLine in donorReadList):
					addedCount = addedCount + 1
				else:
					print("Read " +str(fileLine)+ " was not added.")
				next(gzFile)
				next(gzFile)
				next(gzFile)
	print("Added " +str(addedCount)+ " of " +str(len(donorReadList))+ " variant context donor reads")
	return addedCount

#donorReadList = readDonorReadList(sys.argv[1])
#addedDonorReads = checkDonorReadsAdded(sys.argv[2], donorReadList)
#print("Added " +str(addedDonorReads)+ " of " +str(len(donorReadList))+ " donor reads.")
