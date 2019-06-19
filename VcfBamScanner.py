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
            if sampleid in self.bam_sample_map:
                vcf_to_bam_map[self.vcf_sample_map[sampleid]] = self.bam_sample_map[sampleid]
        return vcf_to_bam_map

    # ===METHODS TO SCAN VCF/BAM FOLDERS AND FILES=============================
    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    def scan_vcf_folders(self, vcffolders):
        self.vaselogger.info("Start scannig VCF files")

        for vcffolder in vcffolders:
            if not os.path.isdir(vcffolder):
                self.vaselogger.info(f"Folder {vcffolder} does not exist.")
                continue
            self.vaselogger.info(f"Scanning VCF files in {vcffolder}")
            for vcf_filename in os.listdir(vcffolder):
                if vcf_filename.endswith((".vcf.gz", ".bcf")):
                    # Do something to scan the file for the sample.
                    try:
                        vcffile = pysam.VariantFile(
                                f"{vcffolder}/{vcf_filename}",
                                "r"
                                )
                        self.vaselogger.debug(
                                "Scanning VCF file "
                                f"{vcffolder}/{vcf_filename}"
                                )
                        if len(vcffile.header.samples) > 0:
                            # This is the sample identifier.
                            sampleid = vcffile.header.samples[0]
                            self.vcf_sample_map[sampleid] = (
                                    f"{vcffolder}/{vcf_filename}"
                                    )
                        vcffile.close()
                    except IOError as ioe:
                        self.vaselogger.warning(
                                f"VCF file {vcffolder}/{vcf_filename} "
                                "does not exist"
                                )
        self.vaselogger.info("Finished scanning VCF files")
        return self.vcf_sample_map

    # Scans the folders containing BAM files and returns a map that
    # links sample ids with bam files.
    def scan_bam_folders(self, bamfolders):
        self.vaselogger.info("Start scanning BAM files")

        # Scan BAM files in all provided folders.
        for bamfolder in bamfolders:
            if not os.path.isdir(bamfolder):
                self.vaselogger.info(f"Folder {bamfolder} does not exist.")
                continue
            self.vaselogger.info(f"Scanning BAM files in {bamfolder}")
            for bamfilename in os.listdir(bamfolder):
                if bamfilename.endswith(".bam"):
                    try:
                        bamfile = pysam.AlignmentFile(
                                f"{bamfolder}/{bamfilename}",
                                "rb"
                                )
                        self.vaselogger.debug(
                                "Scanning BAM file "
                                f"{bamfolder}/{bamfilename}"
                                )
                        if self.bam_has_sample_name(bamfile):
                            # The sample identifier.
                            sampleid = bamfile.header["RG"][0]["SM"]
                            self.bam_sample_map[sampleid] = (
                                    f"{bamfolder}/{bamfilename}"
                                    )
                        bamfile.close()
                    except IOError as ioe:
                        self.vaselogger.warning(
                                f"BAM file {bamfolder}/{bamfilename} "
                                "does not exist"
                                )
        self.vaselogger.info("Finished scanning BAM files")
        return self.bam_sample_map

    # Checks whether the BAM file contains a sample name.
    def bam_has_sample_name(self, bamfile):
        if "RG" in bamfile.header:
            if len(bamfile.header["RG"]) > 0:
                if "SM" in bamfile.header["RG"][0]:
                    return True
        return False
