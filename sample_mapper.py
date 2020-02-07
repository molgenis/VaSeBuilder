"""Builds a list of sample objects.

This module defines Sample objects, consisting of sample id, sample VCF,
sample BAM, and an Argon2 hash of the sample ID. It also defines a method
class, SampleMapper, to build these Sample objects by reading lists of BAM
and VCF files and reading the sample IDs of each. Only samples whose IDs are
found in both a VCF file and a BAM file are stored as Sample objects.
"""

import logging
import os
import subprocess
import sys
import re
import argon2
import pysam


class Sample:
    """Data object that stores metadata about sample file paths.

    Parameters
    ----------
    ID: str
        Sample ID, as found in its SAM 'SM' field or VCF sample column.
    BAM: str
        Path to sample's BAM/CRAM file.
    VCF: str
        Path to sample's VCF/BCF file.
    Hash: str
        Argon2 encoding hash of `ID`.

    Attributes
    ----------
    Hash_ID: str
        Subset from `Hash` containing only the ID hash.
    """

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
        """Return string representation.

        Returns
        -------
        rep_string : str
        """
        rep_string = []
        for key in self.__dict__:
            rep_string.append(f"{key}=\'{self.__dict__[key]}\'")
        rep_string = f"Sample({', '.join(rep_string)})"
        return rep_string

    def get_short_id(self):
        """Return short DNA ID from full ID.

        Returns
        -------
        short_id : str
            Short 6-digit sample ID.
        """
        short_id = re.findall(r"DNA[0-9]{6}", self.ID)[0]
        return short_id


class SampleMapper:
    """Method object to create Sample objects from lists of VCF and BAM files.

    Methods
    -------
    build_sample_maps(bam_list_file, vcf_list_file, make_hash=True)
    """

    def __init__(self):
        return

    @staticmethod
    def read_donor_list_file(list_file):
        """Read in a file containg a list of file paths.

        Parameters
        ----------
        list_file : str

        Returns
        -------
        list
            List of file paths.
        """
        try:
            with open(list_file, "r") as list_file_in:
                file_lines = list_file_in.readlines()
        except IOError as ioe:
            return ioe
        file_lines = list(set(file_lines))
        return [line.strip() for line in file_lines
                if not line.startswith("#")]

    @staticmethod
    def check_donor_files(file_list, file_type):
        """Check existence and valid file types of files in a list of paths.

        Parameters
        ----------
        file_list : list
            List of file paths to check.
        file_type : {'a', 'v'}
            `a` for alignment files, `v` for variant files.

        Returns
        -------
        checked_file_list : list
            List of only file paths leading to existing, valid filetypes.
        """
        if file_type == "a":
            file_types = ("BAM", "CRAM")
        elif file_type == "v":
            file_types = ("VCF", "BCF")
        checked_file_list = []
        warning_list = []
        for file in file_list:
            if not os.path.isfile(file):
                warning_list.append(file)
                continue
            type_check = subprocess.run(
                ["file", "-b", "-z", file],
                stdout=subprocess.PIPE, check=True).stdout.decode()
            if not (file_types[0] in type_check
                    or file_types[1] in type_check):
                warning_list.append(file)
                continue
            checked_file_list.append(file)
        return checked_file_list, warning_list

    @staticmethod
    def read_alignment_sample_ids(bam):
        """Extract sample IDs from a BAM/CRAM file.

        Parameters
        ----------
        bam : str
            Path to BAM/CRAM file.

        Returns
        -------
        list
            List of sample IDs.
        """
        try:
            bamfile = pysam.AlignmentFile(bam, reference_filename=None)
        except IOError as ioe:
            return ioe
        id_list = []
        for rg_line in bamfile.header["RG"]:
            id_list.append(rg_line["SM"])
        return list(set(id_list))

    @staticmethod
    def read_vcf_sample_ids(vcf):
        """Extract sample IDs from a VCF.gz/BCF file.

        Parameters
        ----------
        vcf : str
            Path to VCF.gz/BCF file.

        Returns
        -------
        list
            List of sample IDs.
        """
        try:
            vcffile = pysam.VariantFile(vcf)
        except IOError as ioe:
            return ioe
        return list(set(vcffile.header.samples))

    @staticmethod
    def hash_sample_id(hasher: argon2.PasswordHasher, sample: Sample,
                       remove_id=False):
        """Produce hash of sample.ID.

        Argon2 hash encoding is stored in sample.Hash, and the sample.ID hash
        itself is stored in sample.Hash_ID.

        Parameters
        ----------
        hasher : argon2.PasswordHasher
            Instance of `argon2.PasswordHasher`.
        sample : Sample
            Sample object.
        removeID : bool, optional
            Overwrite sample.ID with sample.Hash_ID if True. Default is False.

        Returns
        -------
        None.
        """
        sample.Hash = hasher.hash(sample.ID)
        sample.Hash_ID = sample.Hash.split("$")[-1]
        if remove_id:
            sample.ID = sample.Hash_ID

    @staticmethod
    def get_alignment_sequence_names(alignment_file):
        """Read the chromosome names from an alignment file.

        Chromosome names are return from a specified BAM or CRAM file by
        checking the header for 'SQ' tags. If present, the SQ entry is checked
        for the 'SN' field, which contains the specific chromosome names.

        Parameters
        ----------
        alignment_file : pysam.AlignmentFile
            Already opened pysam AlignmentFile

        Returns
        -------
        sequence_names : list
            List of chromosome names
        """
        sequence_names = set()
        if "SQ" in alignment_file.header:
            for sn_entry in alignment_file.header["SQ"]:
                sequence_names.add(sn_entry["SN"])
        return sequence_names

    @classmethod
    def build_sample_maps(cls, bam_list_file, vcf_list_file, make_hash=True):
        """Check file lists and produce complete `Sample` objects.

        Wraps other class methods to read in alignment and variant file lists.
        Then checks file paths and extracts sample IDs from them. Then produces
        `Sample` objects from IDs that exist in both alignment and variant
        files. Finally, hashes sample ID using Argon2 if option set to True.

        Parameters
        ----------
        bam_list_file : str
            Path to file containing list of BAM/CRAM files.
        vcf_list_file : str
            Path to file containing list of VCF.gz/BCF files.
        make_hash : bool, optional
            Produces sample.Hash from sample.ID if True. The default is True.

        Returns
        -------
        sample_list : list
            List of `Sample` objects.

        """
        vaselogger = logging.getLogger("VaSe_Logger")

        types = ["a", "v"]
        bam_map = {}
        vcf_map = {}

        for filetype, list_file in zip(types, [bam_list_file, vcf_list_file]):
            files = cls.read_donor_list_file(list_file)
            if isinstance(files, Exception):
                vaselogger.critical("Could not open donor list file:\n"
                                    f"{files}")
                sys.exit()
            files, warnings = cls.check_donor_files(files, filetype)
            if warnings:
                vaselogger.warning("Files do not exist or are not supported "
                                   "formats:\n{}".format("\n".join(warnings)))

            for file in files:
                if filetype == "a":
                    file_ids = cls.read_alignment_sample_ids(file)
                elif filetype == "v":
                    file_ids = cls.read_vcf_sample_ids(file)
                if isinstance(file_ids, Exception):
                    vaselogger.warning(f"Unable to read file: {file_ids}")
                for file_id in file_ids:
                    if filetype == "a":
                        bam_map[file_id] = file
                    if filetype == "v":
                        vcf_map[file_id] = file

# =============================================================================
#         bamlist = cls.read_donor_list_file(bam_list_file)
#         if isinstance(bamlist, Exception):
#             vaselogger.critical("Could not open donor list file:\n"
#                                 f"{bamlist}")
#             sys.exit()
#         bamlist, warnlist = cls.check_donor_files(bamlist, "a")
#         if warnlist:
#             vaselogger.warning("Files do not exist or are not supported "
#                                "alignment files:\n"
#                                "{}".format("\n".join(warnlist)))
#         bam_sample_map = {}
#         for bam in bamlist:
#             bam_ids = cls.read_alignment_sample_ids(bam)
#             if isinstance(bam_ids, Exception):
#                 vaselogger.warning(f"Unable to read file: {bam_ids}")
#             for bam_id in bam_ids:
#                 bam_sample_map[bam_id] = bam
#
#         vcflist = cls.read_donor_list_file(vcf_list_file)
#         if isinstance(vcflist, Exception):
#             vaselogger.critical("Could not open donor list file:\n"
#                                 f"{vcflist}")
#             sys.exit()
#         vcflist, warnlist = cls.check_donor_files(vcflist, "v")
#         if warnlist:
#             vaselogger.warning("Files do not exist or are not supported "
#                                "variant files:\n"
#                                "{}".format("\n".join(warnlist)))
#         vcf_sample_map = {}
#         for vcf in vcflist:
#             vcf_ids = cls.read_vcf_sample_ids(vcf)
#             if isinstance(vcf_ids, Exception):
#                 vaselogger.warning(f"Unable to read file: {vcf_ids}")
#             for vcf_id in vcf_ids:
#                 vcf_sample_map[vcf_id] = vcf
# =============================================================================

        sample_list = []
        for sample_id in list(set(bam_map.keys()) & set(vcf_map.keys())):
            sample_list.append(
                Sample(sample_id, bam_map[sample_id], vcf_map[sample_id])
                )

        if make_hash:
            hasher = argon2.PasswordHasher()
            for sample in sample_list:
                cls.hash_sample_id(hasher, sample)
        elif not make_hash:
            for sample in sample_list:
                sample.Hash_ID = sample.ID
        return sample_list


if __name__ == "__main__":
    SAMPLES = SampleMapper.build_sample_maps(sys.argv[1], sys.argv[2])
