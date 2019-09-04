#!/usr/bin/env python
import logging
import os
import subprocess
import pysam


class VcfBamScanner:
    """VcfBamScanner offers functionality to scan variant (VCF/BCF) and alignment (BAM/CRAM) files.

    The VcfBamScanner can be used to scan folders for variant and alignments files, check whether variant and alignment
    files have sample names and extract them. A list file containing paths to either variant or alignment files can
    also be.scanned. Files within the list files will be a saved per sample name in an identifier."""
    # Constructor that creates two empty hashmaps (dicts).
    def __init__(self):
        """Obtains the vase logger and creates two empty dictionaries to save the variant and alignment files."""
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_sample_map = {}
        self.bam_sample_map = {}
        self.valid_file_types = ["VCF", "BCF", "BAM", "CRAM"]

    # ===METHODS TO GET SAVED DATA FROM THE VCFBAMSCANNER======================
    # Returns the map that links VCF files and samples.
    def get_vcf_sample_map(self):
        """Returns a dictionary with VCF/BCF files saved per sample."""
        return self.vcf_sample_map

    # Returns the map that links BAM files and samples.
    def get_bam_sample_map(self):
        """Returns a dictionary with BAM/CRAM files saved per sample."""
        return self.bam_sample_map

    # Returns a map that links VCF files to BAM files.
    def get_vcf_to_bam_map(self):
        """Returns a dictionary with a variant file linked to an alignment file.

        For each sample, the associated variant file is linked to the associated.alignment file are linked via a
        dictionary. The variant and alignment file are saved as key and value respectively."""
        vcf_to_bam_map = {}
        for sampleid in self.vcf_sample_map:
            if sampleid in self.bam_sample_map:
                vcf_to_bam_map[self.vcf_sample_map[sampleid]] = self.bam_sample_map[sampleid]
        return vcf_to_bam_map

    # ===METHODS TO SCAN VCF/BAM FOLDERS AND FILES=============================
    # Scans the folders containing VCF files and returns a map that
    # links sample ids with vcf files.
    def scan_vcf_folders(self, vcffolders):
        """Scans a list of folders containing variant files and returns a dictionary with valid files per sample name.

        Iterates over a list of paths to folders containing VCF/BCF files. A non-existing folder is skipped. For each
        variant file in a valid folder, the sample name is extracted.and the file is saved under its sample name in a
        dictionary."""
        self.vaselogger.info("Start scannig VCF files")

        for vcffolder in vcffolders:
            if not os.path.isdir(vcffolder):
                self.vaselogger.info(f"Folder {vcffolder} does not exist.")
                continue
            self.vaselogger.info(f"Scanning VCF files in {vcffolder}")
            for vcf_filename in os.listdir(vcffolder):
                if vcf_filename.endswith((".vcf", ".vcf.gz", ".bcf")):
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

    # Scans the folders containing BAM files and returns a map that links sample ids with bam files.
    def scan_bam_folders(self, bamfolders):
        """Scans a list of folders containing alignment files and returns a dictionary with valid files per sample name.

        Iterates over a list of paths to folders containing BAM/CRAM files. Non-existing folders are skipped. For each
        alignment file in a valid folder, the sample name is extracted.and the file is saved under its sample name in a
        dictionary."""
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
        """Reads a list file containing the locations of BAM/CRAM files and saves them per sample."""
        self.bam_sample_map = self.read_donorlistfile(bamcramlistfile)
        return self.bam_sample_map

    # Scan the provided VCF files from a provided donor list file
    def scan_vcf_files(self, vcflistfile):
        """Reads a list file containing the locations of VCF/BCF files and saves them per sample name."""
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
        """Checks and returns whether a BAM/CRAM file has a sample name/identifier."""
        if "RG" in bamfile.header:
            if len(bamfile.header["RG"]) > 0:
                if "SM" in bamfile.header["RG"][0]:
                    return True
        return False

    # Checks whether the VCF file has a sample name
    def vcf_has_sample_name(self, vcffile):
        """Checks and returns whether the VCF/BCF file has a sample name/identifier."""
        if vcffile.header.samples[0] != "":
            return True
        return False

    # Returns the first sample identifier of the VCF file (as for now we only work with one sample per file
    def get_vcf_sample_name(self, vcffileloc):
        """Extracts and returns the list of sample names from a specified VCF file."""
        vcffile = pysam.VariantFile(vcffileloc, "r")
        if self.vcf_has_sample_name(vcffile):
            return list(vcffile.header.samples)
        return None

    # Returns the BAM sample
    def get_bam_sample_name(self, bamfileloc):
        """Extracts the sample name/identifier from a BAM file.

        The sample name is extracted by first checking whether the BAM file contains a sample name via the
        VcfBamScanner method bam_has_sample_name(). If so, the sample name is extracted from the BAM file by returning
        the value of the 'SM' field from the header line starting with '@RG'."""
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
        """Extracts the sample name/identifier from a specified CRAM file.

        The sample name is extracted by first checking whether the CRAM file contains a sample name via the
        VcfBamScanner method cram_has_sample_name(). If so, the sample name is extracted from the CRAM file by
        returning the value of the 'SM' field from the header line starting with '@RG'. None is returned if the CRAM
        file has no sample name."""
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
        """Returns the sample name/identifier for a specified alignemt or variant file."""
        donor_file_type = self.get_donor_file_type(donorfileloc)
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
        """Determines and returns the file type of a specified donor file.

        To determine the file type, the linux command 'file' is used which returns a String containing the file type
        info. The received String is split on spaces and the resulting list is returned."""
        filetype_proc = subprocess.Popen(["file", "-b", "-z", donorfileloc], stdout=subprocess.PIPE)
        filetype_data, filetype_err = filetype_proc.communicate()
        return filetype_data.decode().split(" ")

    # Returns a list of sample identifiers that have both a VCF and a BAM/CRAM
    def get_complete_sample_ids(self):
        """Determines which samples have a variant and alignment file and the list of identifiers."""
        return set(self.vcf_sample_map.keys()) & set(self.bam_sample_map.keys())

    # Extracts the sequence names from an already opened BAM/CRAM file
    def get_alignment_sequence_names(self, alignment_file):
        """Extracts and returns the chromosome names from an alignment file.

        Chromosome names are extracted from a specified BAM or CRAM file by checking the headers for the 'SQ' tag. If
        present the SQ entry is checked for the 'SN' field which contains the specific chromosome names."""
        sequence_names = set()
        if "SQ" in alignment_file.header:
            for sn_entry in alignment_file.header["SQ"]:
                sequence_names.add(sn_entry["SN"])
        return sequence_names

    # Extracts the sequence names from an already opened VCF/BCF file
    def get_variant_sequence_names(self, variant_file):
        """Extracts and returns the chromosome names from a variant file.

        Chromosome names are extracted from a specified VCF or BCF file by obtaining the header contig names via the
        pysam method .header.contigs()."""
        sequence_names = set()
        if len(variant_file.header.contigs) > 0:
            for seqname in list(variant_file.header.contigs):
                sequence_names.add(seqname)
        return sequence_names
