#!/usr/bin/env python
import logging
import os
import subprocess
import pysam


class VcfBamScanner:
    # Constructor that creates two empty hashmaps (dicts).
    def __init__(self):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_sample_map = {}
        self.bam_sample_map = {}
        self.valid_file_types = ["VCF", "BCF", "BAM", "CRAM"]

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

    # Scans the provided BAM/CRAM files from a provided donor list file
    def scan_bamcram_files(self, bamcramlistfile):
        self.bam_sample_map = self.read_donorlistfile(bamcramlistfile)
        return self.bam_sample_map

    # Scan the provided VCF files from a provided donor list file
    def scan_vcf_files(self, vcflistfile):
        self.vcf_sample_map = self.read_donorlistfile(vcflistfile, "v")
        return self.vcf_sample_map

    # Read a specific donor list file
    def read_donorlistfile(self, listfileloc, listtype="a"):
        donor_sample_files = {}
        try:
            with open(listfileloc, "r") as listfile:
                for fileline in listfile:
                    if os.path.isfile(fileline.strip()):
                        sampleid = self.get_sample_id(fileline.strip(), listtype)

                        if sampleid is not None:
                            if type(sampleid) == list:
                                for sid in sampleid:
                                    donor_sample_files[sid] = fileline.strip()
                            else:
                                donor_sample_files[sampleid] = fileline.strip()
        except IOError:
            self.vaselogger.critical(f"Could not open donor list file {listfileloc}")
            exit()
        return donor_sample_files

    # Checks whether the BAM file contains a sample name.
    def bam_has_sample_name(self, bamfile):
        if "RG" in bamfile.header:
            if len(bamfile.header["RG"]) > 0:
                if "SM" in bamfile.header["RG"][0]:
                    return True
        return False

    # Checks whether the VCF file has a sample name
    def vcf_has_sample_name(self, vcffile):
        if vcffile.header.samples[0] != "":
            return True
        return False

    # Returns the first sample identifier of the VCF file (as for now we only work with one sample per file
    def get_vcf_sample_name(self, vcffileloc):
        vcffile = pysam.VariantFile(vcffileloc, "r")
        if self.vcf_has_sample_name(vcffile):
            return list(vcffile.header.samples)
        return None

    # Returns the BAM sample
    def get_bam_sample_name(self, bamfileloc):
        try:
            bamfile = pysam.AlignmentFile(bamfileloc, "rb")
            if self.bam_has_sample_name(bamfile):
                return bamfile.header["RG"][0]["SM"]
            else:
                self.vaselogger.debug(f"BAM file {bamfileloc} has no sample identifier")
                return None
        except IOError:
            self.vaselogger.warning(f"Could not open {bamfileloc} to extract sample identifier")

    # Returns the CRAM sample name
    def get_cram_sample_name(self, cramfileloc):
        try:
            cramfile = pysam.AlignmentFile(cramfileloc, "rc")
            if self.bam_has_sample_name(cramfile):
                return cramfile.header["RG"][0]["SM"]
            else:
                return None
        except IOError:
            self.vaselogger.warning(f"Could not open {cramfileloc} to extract sample identifier")

    # General methods that returns one or more sample ids (in case of a VCF)
    def get_sample_id(self, donorfileloc, donorlisttype):
        donor_file_type = self.get_donor_file_type(donorfileloc)
        print(donor_file_type)
        if donorlisttype == "a":
            if donor_file_type[1] == "BAM":
                return self.get_bam_sample_name(donorfileloc)
            elif donor_file_type[0] == "CRAM":
                return self.get_cram_sample_name(donorfileloc)
        elif donorlisttype == "v":
            if donor_file_type[3] == "(VCF)" or donor_file_type[3] == "(BCF)":
                return self.get_vcf_sample_name(donorfileloc)

    # Returns the filetype of a donor file
    def get_donor_file_type(self, donorfileloc):
        filetype_proc = subprocess.Popen(["file", "-b", "-z", donorfileloc], stdout=subprocess.PIPE)
        filetype_data, filetype_err = filetype_proc.communicate()
        return filetype_data.decode().split(" ")

    # Returns a list of sample identifiers that have both a VCF and a BAM/CRAM
    def get_complete_sample_ids(self):
        return set(self.vcf_sample_map.keys()) & set(self.bam_sample_map.keys())
