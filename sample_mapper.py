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
import argon2
import pysam


class Sample:
    """Data object that stores metadata about sample file paths.

    Parameters
    ----------
    sample_id: str
        Sample ID, as found in its SAM 'SM' field or VCF sample column.
    bam_path: str
        Path to sample's BAM/CRAM file.
    vcf_path: str
        Path to sample's VCF/BCF file.
    argon2_encoding: str
        Argon2 encoding hash of `ID`.

    Attributes
    ----------
    Hash_ID: str
        Subset from `Hash` containing only the ID hash.
    """

    def __init__(self,
                 sample_id: str = None,
                 bam_path: str = None,
                 vcf_path: str = None,
                 argon2_encoding: str = None):
        self.id = sample_id
        self.bam = bam_path
        self.vcf = vcf_path
        self.hash = argon2_encoding
        self.hash_id = None

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
        """Read in a file containing a list of file paths.

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
        """Produce hash of sample.id.

        Argon2 hash encoding is stored in sample.hash, and the sample.id hash
        itself is stored in sample.hash_id.

        Parameters
        ----------
        hasher : argon2.PasswordHasher
            Instance of `argon2.PasswordHasher`.
        sample : Sample
            Sample object.
        removeID : bool, optional
            Overwrite sample.id with sample.hash_id if True. Default is False.

        Returns
        -------
        None.
        """
        sample.hash = hasher.hash(sample.id)
        sample.hash_id = sample.hash.split("$")[-1]
        if remove_id:
            sample.id = sample.hash_id

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

    @staticmethod
    def read_hashtable(hashtable):
        """Read a VaSeBuilder-produced sample ID hashtable of Argon2 encodings."""
        try:
            with open(hashtable) as infile:
                hashes = infile.readlines()
        except IOError as ioe:
            return ioe
        hashes = [line.strip().split("\t") for line in hashes
                  if not line.startswith("#")]
        hashes = {line[0]: line[1] for line in hashes}
        return hashes

    @classmethod
    def build_sample_maps(cls, bams, vcfs, make_hash=True, hashtable=None):
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
            Produces sample.hash from sample.id if True. The default is True.

        Returns
        -------
        sample_list : list
            List of `Sample` objects.

        """
        vaselogger = logging.getLogger("VaSe_Logger")

        types = ["a", "v"]
        bam_map = {}
        vcf_map = {}

        for filetype, files in zip(types, [bams, vcfs]):
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

        sample_list = []
        for sample_id in list(set(bam_map.keys()) & set(vcf_map.keys())):
            sample_list.append(
                Sample(sample_id, bam_map[sample_id], vcf_map[sample_id])
                )

        if make_hash:
            if hashtable is not None:
                hashes = cls.read_hashtable(hashtable)
                if isinstance(hashes, Exception):
                    vaselogger.critical("Unable to read hashtable: {hashes}")
                    sys.exit()
                for sample in sample_list:
                    if sample.id not in hashes:
                        vaselogger.critical("Sample {sample.id} not in hashtable.")
                        sys.exit()
                    sample.hash = hashes[sample.id]
                    sample.hash_id = sample.hash.split("$")[-1]
            elif hashtable is None:
                hasher = argon2.PasswordHasher(memory_cost=1024)
                for sample in sample_list:
                    cls.hash_sample_id(hasher, sample)
        elif not make_hash:
            for sample in sample_list:
                sample.hash_id = sample.id
        return sample_list


if __name__ == "__main__":
    SAMPLES = SampleMapper.build_sample_maps(sys.argv[1], sys.argv[2])
