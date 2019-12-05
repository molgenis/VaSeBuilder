#!/usr/bin/env python
# import logging
import os
import subprocess
import pysam
import sys
sys.path.append("/groups/umcg-atd/tmp03/umcg-tmedina/repos/PyPackages/Argon2")
import argon2
# from dataclasses import dataclass


# =============================================================================
# @dataclass
# class Sample:
#     ID: str = "NULL"
#     BAM: str = "NULL"
#     VCF: str = "NULL"
#     Hash: str = "NULL"
# =============================================================================


class Sample:
    def __init__(self,
                 ID: str = "NULL",
                 BAM: str = "NULL",
                 VCF: str = "NULL",
                 Hash: str = "NULL"):
        self.ID = ID
        self.BAM = BAM
        self.VCF = VCF
        self.Hash = Hash

    def __repr__(self):
        rep_string = []
        for key in self.__dict__:
            rep_string.append(f"{key}=\'{self.__dict__[key]}\'")
        rep_string = f"Sample({', '.join(rep_string)})"
        return rep_string


class SampleMapper:
    """VcfBamScanner offers functionality to scan variant (VCF/BCF) and alignment (BAM/CRAM) files.

    The VcfBamScanner can be used to scan folders for variant and alignments files, check whether variant and alignment
    files have sample names and extract them. A list file containing paths to either variant or alignment files can
    also be.scanned. Files within the list files will be a saved per sample name in an identifier.

    Attributes
    ----------
    vcf_sample_map : dict
        Dictionary with variant files per sample
    bam_sample_map : dict
        Dictionary with alignment files per sample
    valid_file_types : list of str
        Valid variant and alignment file types
    """

    def __init__(self):
        """DOCSTRING"""
        return

    @staticmethod
    def read_donor_list_file(list_file):
        try:
            with open(list_file, "r") as list_file_in:
                file_lines = list_file_in.readlines()
        except IOError:
            # self.vaselogger.critical(f"Could not open donor list file: {list_file}")
            print(f"Could not open donor list file: {list_file}")
            exit()
        file_lines = list(set(file_lines))
        return [line.strip() for line in file_lines if not line.startswith("#")]

    @staticmethod
    def check_donor_files(file_list, file_type):
        if file_type == "a":
            file_types = ("BAM", "CRAM")
        elif file_type == "v":
            file_types = ("VCF", "BCF")
        checked_file_list = []
        for file in file_list:
            if not os.path.isfile(file):
                # self.vaselogger.warning(f"Listed path {file} is not a file.")
                print(f"Listed path {file} is not a file.")
                continue
            type_check = subprocess.run(["file", "-b", "-z", file], stdout=subprocess.PIPE).stdout.decode()
            if not (file_types[0] in type_check or file_types[1] in type_check):
                # self.vaselogger.warning(f"Listed path {alignment} is not a {file_types[0]} or {file_types[1]} file.")
                print(f"Listed path {file} is not a {file_types[0]} or {file_types[1]} file.")
                continue
            checked_file_list.append(file)
        return checked_file_list

    @staticmethod
    def read_alignment_sample_ids(bam):
        try:
            bamfile = pysam.AlignmentFile(bam, reference_filename=None)
        except IOError:
            # self.vaselogger.warning(f"Unable to read file: {bam}")
            print(f"Unable to read file: {bam}")
            return
        id_list = []
        for rg_line in bamfile.header["RG"]:
            id_list.append(rg_line["SM"])
        return list(set(id_list))

    @staticmethod
    def read_vcf_sample_ids(vcf):
        try:
            vcffile = pysam.VariantFile(vcf)
        except IOError:
            # self.vaselogger.warning(f"Unable to read file: {vcf}")
            print(f"Unable to read file: {vcf}")
            return
        return list(set(vcffile.header.samples))

    @staticmethod
    def hash_sample_id(hasher: argon2.PasswordHasher, sample: Sample, removeID=False):
        sample.Hash = hasher.hash(sample.ID)
        sample.Hash_ID = sample.Hash.split("$")[-1]
        if removeID:
            sample.ID = sample.Hash_ID
        return

    @classmethod
    def build_sample_maps(cls, bam_list_file, vcf_list_file, make_hash=True):
        bamlist = cls.read_donor_list_file(bam_list_file)
        bamlist = cls.check_donor_files(bamlist, "a")
        bam_sample_map = {}
        for bam in bamlist:
            bam_ids = cls.read_alignment_sample_ids(bam)
            for bam_id in bam_ids:
                bam_sample_map[bam_id] = bam

        vcflist = cls.read_donor_list_file(vcf_list_file)
        vcflist = cls.check_donor_files(vcflist, "v")
        vcf_sample_map = {}
        for vcf in vcflist:
            vcf_ids = cls.read_vcf_sample_ids(vcf)
            for vcf_id in vcf_ids:
                vcf_sample_map[vcf_id] = vcf

        sample_list = []
        for sample_id in list(set(bam_sample_map.keys()) & set(vcf_sample_map.keys())):
            sample_list.append(Sample(sample_id, bam_sample_map[sample_id], vcf_sample_map[sample_id]))

        if make_hash:
            hasher = argon2.PasswordHasher()
            for sample in sample_list:
                cls.hash_sample_id(hasher, sample)
        return sample_list


if __name__ == "__main__":
    samples = SampleMapper.build_sample_maps(sys.argv[1], sys.argv[2])
