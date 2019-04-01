import sys
import gzip

# Reads the list of acceptor reads.
def readAcceptorReadsList(acceptorReadFile):
	acceptorReads = []
	with open(acceptorReadFile, 'r') as arFile:
		next(arFile)	# Skip the header line
		for fileLine in arFile:
			fileLine = fileLine.strip()
			fileLineData = fileLine.split("\t")
			acceptorReads.extend(fileLineData[1:])
	return acceptorReads


# Checks whether the identified acceptor reads have been removed.
def checkAcceptorReadsRemoved(gzResultsFile, acceptorReadFile):
	acceptorReadList = readAcceptorReadsList()
	removalCount = len(acceptorReadList)
	
	# Check if all accepor reads have been removed from the VaSe fastq file.
	with gzip.open(gzResultsFile, 'rt') as gzFile:
		for fileLine in gzFile:
			fileLine = fileLine.strip()
			if(fileLine.startswith('@')):
				if(fileLine[1:-2] in acceptorReadList):
					print("Read " +str(fileLine[1:-2])+ " was not removed")
					removalCount = removalCount - 1
	return removalCount

#acceptorReadList = readAcceptorReadsList(sys.argv[1])
#removedAcceptorReads = checkAcceptorReadsRemoved(sys.argv[2], acceptorReadList)
#print("Removed " +str(removedAcceptorReads)+ " of " +str(len(acceptorReadList))+ " reads")
