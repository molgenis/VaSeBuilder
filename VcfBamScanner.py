#!/usr/bin/env python
import logging
import os
import pysam


class VcfBamScanner:
    # Constructor that creates two empty hashmaps (dicts).
    def __init__(self):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_sample_map = {}
        self.bam_sample_map = {}

    # ===METHODS TO GET SAVED DATA FROM THE VCFBAMSCANNER======================
    # Returns the map that links VCF files and samples.
    def get_vcf_sample_map(self):
        return self.vcf_sample_map

    # Returns the map that links BAM files and samples.
    def get_bam_sample_map(self):
        return self.bam_sample_map

    # Returns a map that links VCF files to BAM files.
    def get_vcf_to_bam_map(self):
        vcf_to_bam_map = {}
        for sampleid in self.vcf_sample_map:
            if (sampleid in self.bam_sample_map):
                vcf_to_bam_map[self.vcf_sample_map[sampleid]] = self.bam_sample_map[sampleid]
        return vcf_to_bam_map

    # ===METHODS TO SCAN VCF/BAM FOLDERS AND FILES=============================
    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    def scan_vcf_folders(self, vcffolders):
        self.vaselogger.info("Start scannig VCF files")

        for vcffolder in vcffolders:
            if (os.path.isdir(vcffolder)):
                self.vaselogger.info("Scanning VCF files in " + vcffolder)
                for vcf_filename in os.listdir(vcffolder):
                    if (vcf_filename.endswith(".vcf")
                       or vcf_filename.endswith(".vcf.gz")):
                        # Do something to scan the file for the sample.
                        try:
                            vcfFile = pysam.VariantFile(
                                    vcffolder + "/" + vcf_filename,
                                    'r'
                                    )
                            self.vaselogger.debug(
                                    "Scanning VCF file "
                                    + vcffolder + "/" + vcf_filename
                                    )
                            if (len(vcfFile.header.samples) > 0):
                                # This is the sample identifier.
                                sampleid = vcfFile.header.samples[0]
                                self.vcf_sample_map[sampleid] = (vcffolder
                                                                 + "/"
                                                                 + vcf_filename)
                            vcfFile.close()
                        except IOError as ioe:
                            self.vaselogger.warning(
                                    "VCF file "
                                    + vcffolder + "/" + vcf_filename
                                    + " does not exist"
                                    )
            else:
                self.vaselogger.info(
                        "Folder " + vcffolder + " does not exist."
                        )
        self.vaselogger.info("Finished scanning VCF files")
        return self.vcf_sample_map

    # Scans the folders containing BAM files and returns a map that
    # links sample ids with bam files.
    def scan_bam_folders(self, bamFolders):
        self.vaselogger.info("Start scanning BAM files")

        # Scan BAM files in all provided folders.
        for bamFolder in bamFolders:
            if (os.path.isdir(bamFolder)):
                self.vaselogger.info("Scanning BAM files in " + bamFolder)
                for bamFileName in os.listdir(bamFolder):
                    if (bamFileName.endswith(".bam")):
                        try:
                            bamFile = pysam.AlignmentFile(
                                    bamFolder + "/" + bamFileName,
                                    'rb'
                                    )
                            self.vaselogger.debug(
                                    "Scanning BAM file "
                                    + bamFolder + "/" + bamFileName
                                    )
                            if (self.bam_has_sample_name(bamFile)):
                                # The sample identifier.
                                sampleid = bamFile.header["RG"][0]["SM"]
                                self.bam_sample_map[sampleid] = (
                                        bamFolder + "/" + bamFileName
                                        )
                            bamFile.close()
                        except IOError as ioe:
                            self.vaselogger.warning(
                                    "BAM file "
                                    + bamFolder + "/" + bamFileName
                                    + " does not exist"
                                    )
            else:
                self.vaselogger.info(
                        "Folder " + bamFolder + " does not exist."
                        )
        self.vaselogger.info("Finished scanning BAM files")
        return self.bam_sample_map

    # Checks whether the BAM file contains a sample name.
    def bam_has_sample_name(self, bamFile):
        if ("RG" in bamFile.header):
            if (len(bamFile.header["RG"]) > 0):
                if ("SM" in bamFile.header["RG"][0]):
                    return True
        return False
