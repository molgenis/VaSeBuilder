#!/usr/bin/env python
import logging
import os
import pysam


class VcfBamScanner:
    # Constructor that creates two empty hashmaps (dicts).
    def __init__(self):
        self.vaseLogger = logging.getLogger("VaSe_Logger")
        self.vcfSampleMap = {}
        self.bamSampleMap = {}

    # ===METHODS TO GET SAVED DATA FROM THE VCFBAMSCANNER======================
    # Returns the map that links VCF files and samples.
    def getVcfSampleMap(self):
        return self.vcfSampleMap

    # Returns the map that links BAM files and samples.
    def getBamSampleMap(self):
        return self.bamSampleMap

    # Returns a map that links VCF files to BAM files.
    def getVcfToBamMap(self):
        vcfToBamMap = {}
        for sampleid in self.vcfSampleMap:
            if (sampleid in self.bamSampleMap):
                vcfToBamMap[self.vcfSampleMap[sampleid]] = self.bamSampleMap[sampleid]
        return vcfToBamMap

    # ===METHODS TO SCAN VCF/BAM FOLDERS AND FILES=============================
    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    def scanVcfFolders(self, vcfFolders):
        self.vaseLogger.info("Start scannig VCF files")

        for vcfFolder in vcfFolders:
            if (os.path.isdir(vcfFolder)):
                self.vaseLogger.info(f"Scanning VCF files in {vcfFolder}")
                for vcfFileName in os.listdir(vcfFolder):
                    if (vcfFileName.endswith(".vcf")
                       or vcfFileName.endswith(".vcf.gz")):
                        # Do something to scan the file for the sample.
                        try:
                            vcfFile = pysam.VariantFile(
                                    f"{vcfFolder}/{vcfFileName}",
                                    "r"
                                    )
                            self.vaseLogger.debug(
                                    "Scanning VCF file "
                                    f"{vcfFolder}/{vcfFileName}"
                                    )
                            if (len(vcfFile.header.samples) > 0):
                                # This is the sample identifier.
                                sampleid = vcfFile.header.samples[0]
                                self.vcfSampleMap[sampleid] = (
                                        f"{vcfFolder}/{vcfFileName}"
                                        )
                            vcfFile.close()
                        except IOError as ioe:
                            self.vaseLogger.warning(
                                    f"VCF file {vcfFolder}/{vcfFileName} "
                                    "does not exist"
                                    )
            else:
                self.vaseLogger.info(f"Folder {vcfFolder} does not exist.")
        self.vaseLogger.info("Finished scanning VCF files")
        return self.vcfSampleMap

    # Scans the folders containing BAM files and returns a map that
    # links sample ids with bam files.
    def scanBamFolders(self, bamFolders):
        self.vaseLogger.info("Start scanning BAM files")

        # Scan BAM files in all provided folders.
        for bamFolder in bamFolders:
            if (os.path.isdir(bamFolder)):
                self.vaseLogger.info(f"Scanning BAM files in {bamFolder}")
                for bamFileName in os.listdir(bamFolder):
                    if (bamFileName.endswith(".bam")):
                        try:
                            bamFile = pysam.AlignmentFile(
                                    f"{bamFolder}/{bamFileName}",
                                    "rb"
                                    )
                            self.vaseLogger.debug(
                                    "Scanning BAM file "
                                    f"{bamFolder}/{bamFileName}"
                                    )
                            if (self.bamHasSampleName(bamFile)):
                                # The sample identifier.
                                sampleid = bamFile.header["RG"][0]["SM"]
                                self.bamSampleMap[sampleid] = (
                                        f"{bamFolder}/{bamFileName}"
                                        )
                            bamFile.close()
                        except IOError as ioe:
                            self.vaseLogger.warning(
                                    f"BAM file {bamFolder}/{bamFileName} "
                                    "does not exist"
                                    )
            else:
                self.vaseLogger.info(f"Folder {bamFolder} does not exist.")

        self.vaseLogger.info("Finished scanning BAM files")
        return self.bamSampleMap

    # Checks whether the BAM file contains a sample name.
    def bamHasSampleName(self, bamFile):
        if ("RG" in bamFile.header):
            if (len(bamFile.header["RG"]) > 0):
                if ("SM" in bamFile.header["RG"][0]):
                    return True
        return False
