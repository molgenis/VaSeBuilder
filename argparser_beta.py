# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 21:12:05 2020

@author: tdmedina
"""

import argparse
import os
import sys
import subprocess


class MyArgumentParser(argparse.ArgumentParser):
    def convert_arg_line_to_args(self, arg_line):
        return arg_line.split()


class CustomHelp(argparse.HelpFormatter):
    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
        else:
            parts = []
            if action.nargs == 0:
                parts.extend(action.option_strings)
            else:
                default = self._get_default_metavar_for_optional(action)
                args_string = self._format_args(action, default)
                for option_string in action.option_strings:
                    parts.append("%s" % (option_string))
                parts[-1] += " %s " % args_string
            return ', '.join(parts)


class VctorArgs:
    def __init__(self):
        self.vctor_parser = self.make_vctor_parser()

    def make_vctor_parser(cls):
        vctor_args = MyArgumentParser(formatter_class=CustomHelp)
        subparsers = vctor_args.add_subparsers(
            dest="runmode", required=True,
            metavar="{BuildSpikeIns, AssembleValidationSet, BuildValidationSet}"
            )

        # ===Equivalent to D, DC, P, PC, and X modes================================================
        parser_spike = subparsers.add_parser(name="BuildSpikeIns", formatter_class=CustomHelp)
        parser_spike.add_argument("-m", "--output-mode", required=True, choices=["A", "D", "P"],
                                  help=("How to produce outputs. "
                                        "A: Output one VCF, BAM, and variant context file with all variant contexts; "
                                        "D: Output VCF, BAM, and variant context files per sample.; "
                                        "P: Output VCF, BAM, and variant context files per variant context."))
        # Make varcons using acceptor BAM or use existing varcon file(s).
        template_arg = parser_spike.add_mutually_exclusive_group(required=True)
        template_arg.add_argument("-c", "--varcon", nargs="+", dest="varcon_in",
                                  type=cls.is_existing_file, metavar=("<varconfile>", "<varconfile2>"),
                                  help="Pre-made variant context file(s).")
        template_arg.add_argument("-cL", "--varcon-list", dest="varcon_in",
                                  type=cls.are_existing_files, metavar="<file>",
                                  help="Pre-made variant context files listed per line in <file>.")
        template_arg.add_argument("-a", "--acceptor-bam",
                                  type=cls.is_alignment_file, metavar="<bam>",
                                  help="Acceptor BAM or CRAM file.")
        # Donor BAM file(s).
        bam_arg = parser_spike.add_mutually_exclusive_group(required=True)
        bam_arg.add_argument("-b", "--donor-bam", nargs="+", dest="donor_bam",
                             type=cls.is_alignment_file, metavar=("<bam>", "<bam2>"),
                             help="Donor BAM or CRAM file(s).")
        bam_arg.add_argument("-bL", "--donor-bam-list", dest="donor_bam",
                             type=cls.are_alignment_files, metavar="<file>",
                             help="Donor BAM or CRAM files listed per line in <file>.")
        # Donor VCF file(s).
        vcf_arg = parser_spike.add_mutually_exclusive_group(required=True)
        vcf_arg.add_argument("-v", "--donor-vcf", nargs="+", dest="donor_vcf",
                             type=cls.is_variant_file, metavar=("<vcf>", "<vcf2>"),
                             help="Donor VCF file(s).")
        vcf_arg.add_argument("-vL", "--donor-vcf-list", dest="donor_vcf",
                             type=cls.are_variant_files, metavar="<file>",
                             help="Donor VCF files listed per line in <file>.")
        # Optionals.
        parser_spike.add_argument("--no-hash", action="store_true",
                                  help="Use original sample IDs without hashing with Argon2.")
        parser_spike.add_argument("--varcon-only", action="store_true",
                                  help="Suppress BAM and VCF output and only output variant context file(s).")
        parser_spike.add_argument("--no-merge", action="store_true",
                                  help="Do not merge overlapping contexts from the same sample.")
        parser_spike.add_argument("--suppress-conflict-check", action="store_true",
                                  help="Ignore conflicts between contexts from different samples (UNSTABLE).")

        # ===Parent parser for the two variant set building parsers=================================
        validation_parent = subparsers.add_parser(name="_parentvalset",
                                                  formatter_class=CustomHelp,
                                                  add_help=False)
        # Acceptor FastQ args.
        fq1_arg = validation_parent.add_mutually_exclusive_group(required=True)
        fq1_arg.add_argument("-1", "--acceptor-fq-r1", nargs="+", dest="acceptor_fq_1",
                             type=cls.is_existing_file, metavar=("<fastqR1>", "<fastqR1_2>"),
                             help="Acceptor FastQ R1 file(s)")
        fq1_arg.add_argument("-1L", "--acceptor-fq-r1-list", dest="acceptor_fq_1",
                             type=cls.are_existing_files, metavar="<file>",
                             help="Acceptor FastQ R1 files listed per line in <file>.")
        fq2_arg = validation_parent.add_mutually_exclusive_group(required=True)
        fq2_arg.add_argument("-2", "--acceptor-fq-r2", nargs="+", dest="acceptor_fq_2",
                             type=cls.is_existing_file, metavar=("<fastqR2>", "<fastqR2_2>"),
                             help="Acceptor FastQ R2 file(s)")
        fq2_arg.add_argument("-2L", "--acceptor-fq-r2-list", dest="acceptor_fq_2",
                             type=cls.are_existing_files, metavar="<file>",
                             help="Acceptor FastQ R2 files listed per line in <file>.")
        # Optionals.
        validation_parent.add_argument("--seed", default=2,
                                       type=int, metavar="<int>",
                                       help="Random seed used to randomly distribute spike-in reads.")
        validation_parent.add_argument("-av", "--acceptor-vcf",
                                       type=cls.is_variant_file, metavar="<vcf>",
                                       help="Acceptor VCF file, used to make hybrid validation VCF.")

        # ===Equivalent to AC and AB modes==========================================================
        parser_assemble = subparsers.add_parser(name="AssembleValidationSet",
                                                formatter_class=CustomHelp,
                                                parents=[validation_parent])
        # Varcons xor varcon list.
        vacon_arg = parser_assemble.add_mutually_exclusive_group(required=True)
        vacon_arg.add_argument("-c", "--varcon", nargs="+", dest="varcon_in",
                               type=cls.is_existing_file, metavar=("<varconfile>", "<varconfile2>"),
                               help="Pre-made variant context file(s).")
        vacon_arg.add_argument("-cL", "--varcon-list", dest="varcon_in",
                               type=cls.are_existing_files, metavar="<file>",
                               help="Pre-made variant context files listed per line in <file>.")
        # Spike-in read files.
        spike_read_args = parser_assemble.add_mutually_exclusive_group(required=True)
        spike_read_args.add_argument("-kb", "--spike-in-bam", nargs="+", dest="spike_in_bam",
                                     type=cls.is_alignment_file, metavar=("<bam>", "<bam2>"),
                                     help="Pre-built spike-in BAM file(s).")
        spike_read_args.add_argument("-kbL", "--spike-in-bam-list", dest="spike_in_bam",
                                     type=cls.are_alignment_files, metavar="<file>",
                                     help="Pre-built spike-in BAM files listed per line in <file>.")
        spike_read_args.add_argument("-kfq", "--spike-in-fastq-list", dest="spike_in_fastq",
                                     type=cls.are_existing_fastqs, metavar="<file>",
                                     help="Pre-built spike-in FastQ files with pairs listed tab-separated per line in <file>.")
        # Spike-in VCF files (optional).
        spike_vcf_args = parser_assemble.add_mutually_exclusive_group()
        spike_vcf_args.add_argument("-kv", "--spike-in-vcf", nargs="+", dest="spike_in_vcf",
                                    type=cls.is_variant_file, metavar=("<vcf>", "<vcf2>"),
                                    help="Pre-built spike-in VCF file(s).")
        spike_vcf_args.add_argument("-kvL", "--spike-in-vcf-list", dest="spike_in_vcf",
                                    type=cls.are_variant_files, metavar="<file>",
                                    help="Pre-built spike-in VCF files listed per line in <file>.")

        # ===Equivalent to F mode===================================================================
        parser_full = subparsers.add_parser(name="BuildValidationSet", formatter_class=CustomHelp,
                                            parents=[validation_parent])
        parser_full.add_argument("-a", "--acceptor-bam", required=True,
                                 type=cls.is_alignment_file, metavar="<bam>",
                                 help="Acceptor BAM or CRAM file.")
        # Donor BAM file(s).
        bam_arg = parser_full.add_mutually_exclusive_group(required=True)
        bam_arg.add_argument("-b", "--donor-bam", nargs="+", dest="donor_bam",
                             type=cls.is_alignment_file, metavar=("<bam>", "<bam2>"),
                             help="Donor BAM or CRAM file(s).")
        bam_arg.add_argument("-bL", "--donor-bam-list", dest="donor_bam",
                             type=cls.are_alignment_files, metavar="<file>",
                             help="Donor BAM or CRAM files listed per line in <file>.")
        # Donor VCF file(s).
        vcf_arg = parser_full.add_mutually_exclusive_group(required=True)
        vcf_arg.add_argument("-v", "--donor-vcf", nargs="+", dest="donor_vcf",
                             type=cls.is_variant_file, metavar=("<vcf>", "<vcf2>"),
                             help="Donor VCF file(s).")
        vcf_arg.add_argument("-vL", "--donor-vcf-list", dest="donor_vcf",
                             type=cls.are_variant_files, metavar="<file>",
                             help="Donor VCF files listed per line in <file>.")
        # Optionals.
        parser_full.add_argument("--no-hash", action="store_true",
                                 help="Use original sample IDs without hashing with Argon2.")

        # ===Universal options======================================================================
        vctor_args.add_argument("-V", "--version", action="version", version="Vctor v.0.1")
        vctor_args.add_argument("-r", "--reference",
                                type=cls.is_existing_file, metavar="<fasta>",
                                help="Reference sequence fasta")
        vctor_args.add_argument("-o", "--out-dir", default="./",
                                type=cls.is_valid_directory, metavar="<path>",
                                help="Output directory")
        vctor_args.add_argument("-l", "--log", metavar="<str>",
                                help="Log output file name")
        vctor_args.add_argument("--debug", action="store_true",
                                help="Log with maximum verbosity")
        return vctor_args

    @classmethod
    def is_alignment_file(cls, file):
        cls.is_existing_file(file)
        type_check = subprocess.run(
            ["file", "-b", "-z", file],
            stdout=subprocess.PIPE, check=True
            ).stdout.decode()
        if not ("BAM" in type_check or "CRAM" in type_check):
            raise argparse.ArgumentTypeError(f"File {file} is not a supported"
                                             " alignment file type.")
        return file

    @classmethod
    def are_alignment_files(cls, listfile):
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_alignment_file(test_file)
        return files

    @classmethod
    def is_variant_file(cls, file):
        cls.is_existing_file(file)
        type_check = subprocess.run(
            ["file", "-b", "-z", file],
            stdout=subprocess.PIPE, check=True
            ).stdout.decode()
        if not ("VCF" in type_check or "BCF" in type_check):
            raise argparse.ArgumentTypeError(f"File {file} is not a supported"
                                             " variant call file type.")
        return file

    @classmethod
    def are_variant_files(cls, listfile):
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_variant_file(test_file)
        return files

    @staticmethod
    def is_existing_file(file):
        if not os.path.isfile(file):
            raise argparse.ArgumentTypeError(f"File {file} does not exist.")
        return file

    @classmethod
    def are_existing_files(cls, listfile):
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_existing_file(test_file)
        return files

    @classmethod
    def are_existing_fastqs(cls, fastq_pair_file):
        cls.is_existing_file(fastq_pair_file)
        with open(fastq_pair_file) as infile:
            file_pairs = infile.readlines()
        file_pairs = [pair.strip().split("\t") for pair in file_pairs
                      if not pair.startswith("#")]
        for fq1, fq2 in file_pairs:
            cls.is_existing_file(fq1)
            cls.is_existing_file(fq2)
        return file_pairs

    @staticmethod
    def is_valid_directory(directory):
        realpath = os.path.abspath(directory)
        if not (os.path.isdir(realpath) and os.access(realpath, os.X_OK | os.W_OK)):
            raise argparse.ArgumentTypeError(f"Directory {directory} does not exist "
                                             "or is not writeable.")
        return realpath
