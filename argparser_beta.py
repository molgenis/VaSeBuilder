# -*- coding: utf-8 -*-
"""Argument Parser for VaSe

Can be run by initializing VaSeParser class, then running setup().

Created on Tue Feb 18 21:12:05 2020
@author: tdmedina
"""

import argparse
import os
import subprocess
import datetime
import pysam
# from collections import OrderedDict


class CustomHelp(argparse.HelpFormatter):
    """Custom help formatter_class that only displays metavar once."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            default = self._get_default_metavar_for_positional(action)
            metavar, = self._metavar_formatter(action, default)(1)
            return metavar
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

    def _format_action(self, action):
        parts = super()._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts


class VaSeParser(argparse.ArgumentParser):
    """Custom ArgumentParser class with pre-built setup method."""

    def convert_arg_line_to_args(self, arg_line):
        """Read arguments from file.

        This method overwrites the ArgumentParser default. File arguments
        are expected to be --long-name=arg1,arg2,arg3 per line.
        """
        if arg_line.startswith("#"):
            return []
        arg_line = arg_line.strip().split(" ")
        arg_line= [arg for arg in arg_line if arg]
        # arg_line = [y.strip() for x in arg_line for y in x.split(",")]
        return arg_line

    def setup(self):
        """Set up the argument parser for VaSeBuilder."""
        self.formatter_class = CustomHelp
        self.fromfile_prefix_chars = "@"

        subparsers = self.add_subparsers(
            title="Subcommands",
            dest="runmode",  # required=True, <-- This only works for Py3.7+
            metavar="{BuildSpikeIns | AssembleValidationSet | BuildValidationSet}"
            )

        self.add_argument("-V", "--version", action="version", version="VaSe v.0.1")

        # ===Universal options=================================================
        univ_parent = subparsers.add_parser(name="_universalparent",
                                            formatter_class=CustomHelp,
                                            add_help=False)
        universals = univ_parent.add_argument_group(title="Universal Options")
        universals.add_argument("-r", "--reference", required=True,
                                type=self.is_existing_file, metavar="<fasta>",
                                help="Reference sequence fasta")
        universals.add_argument("-o", "--out-dir", default="./",
                                type=self.is_valid_directory, metavar="<path>",
                                help="Output directory")
        universals.add_argument("-l", "--log", metavar="<str>",
                                help="Log output file name")
        universals.add_argument("--debug", action="store_true",
                                help="Log with maximum verbosity")

        # ===Parent parser for context building modes.=========================
        context_parent = subparsers.add_parser(name="_contextparent",
                                               formatter_class=CustomHelp,
                                               add_help=False)
        # Donor BAM file(s).
        bam_arg = context_parent.add_mutually_exclusive_group(required=True)
        bam_arg.add_argument("-b", "--donor-bam", nargs="+", dest="donor_bams",
                             type=self.is_alignment_file, metavar=("<bam>", "<bam2>"),
                             help="Donor BAM or CRAM file(s).")
        bam_arg.add_argument("-bL", "--donor-bam-list", dest="donor_bams",
                             type=self.are_alignment_files, metavar="<file>",
                             help="Donor BAM or CRAM files listed per line in <file>.")
        # Donor VCF file(s).
        vcf_arg = context_parent.add_mutually_exclusive_group(required=True)
        vcf_arg.add_argument("-v", "--donor-vcf", nargs="+", dest="donor_vcfs",
                             type=self.is_variant_file, metavar=("<vcf>", "<vcf2>"),
                             help="Donor VCF file(s).")
        vcf_arg.add_argument("-vL", "--donor-vcf-list", dest="donor_vcfs",
                             type=self.are_variant_files, metavar="<file>",
                             help="Donor VCF files listed per line in <file>.")
        # Optionals relating to filtering.
        filter_controls = context_parent.add_argument_group("Filtering Controls")
        filter_controls.add_argument("-f", "--inclusion-filter",
                                     type=self.is_valid_filter_file, metavar="<file>",
                                     help=("List of variants to include, tab-separated as <sample><chrom><pos><ref><alt>."
                                           "Subsetting this list can be done with -s option."))
        filter_controls.add_argument("-s", "--subset-filter", nargs="+",
                                     type=self.is_valid_filter_format,
                                     metavar=("<Column>:<value>[,<value2>,...]",
                                              "<Column2>:<value>[,<value2>,...]"),
                                     help="Specify subset inclusion filter criteria by column and values.")
        filter_controls.add_argument("-p", "--prioritization", nargs="+",
                                     type=self.is_valid_filter_format,
                                     metavar=("<Column>:<value>[,<value2>,...]",
                                              "<Column2>:<value>[,<value2>,...]"),
                                     help=("In the event of a variant context conflict, prioritize inclusion of "
                                           "variants by values in column. Multiple columns should be added from "
                                           "highest to lowest priority, as should values in those columns. Values "
                                           "not mentioned will be assigned as lowest priority, and non-existant "
                                           "columns will be ignored, but will show a warning. Ex: "
                                           "ImportantColumn:best,medium,worst LessImportant:worst"))
        # Optionals related to context creation.
        context_controls = context_parent.add_argument_group("Context Controls")
        context_controls.add_argument("--no-merge", dest="merge", action="store_false",
                                      help="Do not merge overlapping contexts from the same sample. (FUTURE)")
        context_controls.add_argument("--suppress-conflict-check", action="store_true",
                                      help="Ignore conflicts between contexts from different samples. (FUTURE)")
        context_controls.add_argument("--add-secondary-variants", action="store_true",
                                      help=("When using any kind of variant filtering, if an excluded variant "
                                            "overlaps an included variant context, include it in the VCF output. (FUTURE)"))
        context_controls.add_argument("-vo", "--varcon-out",
                                      default="VaSe_" + str(datetime.date.today()) + ".varcon",
                                      metavar="<str>",
                                      help="Output variant context file name. (Default='VaSe_<date>.varcon')")
        # Options, misc.
        hashing = context_parent.add_mutually_exclusive_group()
        hashing.add_argument("--no-hash", dest="make_hash", action="store_false",
                             help="Use original sample IDs without hashing with Argon2.")
        hashing.add_argument("-x", "--hashtable", type=self.is_existing_file, metavar="<file>",
                             help="Use existing VaSeBuilder-made hashtable file to replace sample IDs.")

        # ===Equivalent to D, DC, P, PC, and X modes================================================
        parser_spike = subparsers.add_parser(name="BuildSpikeIns",
                                             formatter_class=CustomHelp,
                                             parents=[univ_parent, context_parent],
                                             help=("Build variant contexts and spike-ins. "
                                                   "Alternatively, can output variant contexts only, "
                                                   "or build spike-ins from pre-made contexts."))
        parser_spike.add_argument("-m", "--output-mode", required=True, choices=["A", "D", "P", "V"],
                                  help=("How to produce outputs. "
                                        "A: Output one VCF and one BAM file for all variant contexts; "
                                        "D: Output one VCF and BAM file per sample. (FUTURE); "
                                        "P: Output one VCF and BAM file per variant context."
                                        "V: Output variant context file only."))
        # Make varcons using acceptor BAM or use existing varcon file(s).
        template_arg = parser_spike.add_mutually_exclusive_group(required=True)
        template_arg.add_argument("-c", "--varcon", dest="varcons_in",
                                  type=self.is_existing_file, metavar="<varconfile>",
                                  help="Pre-made variant context file.")
        # template_arg.add_argument("-c", "--varcon", nargs="+", dest="varcons_in",
        #                           type=self.is_existing_file, metavar=("<varconfile>", "<varconfile2>"),
        #                           help="Pre-made variant context file(s).")
        # template_arg.add_argument("-cL", "--varcon-list", dest="varcons_in",
        #                           type=self.are_existing_files, metavar="<file>",
        #                           help="Pre-made variant context files listed per line in <file>.")
        template_arg.add_argument("-a", "--acceptor-bam",
                                  type=self.is_alignment_file, metavar="<bam>",
                                  help="Acceptor BAM or CRAM file.")
        # Optionals.
        # output_type.add_argument("--varcon-only", action="store_true",
        #                          help="Suppress BAM and VCF output and only output variant context file(s).")
        parser_spike.add_argument("-O", "--output-type", choices=["B", "U", "F"], default="B",
                                  help="Output <B>am, <U>bam, or <F>astQ spike-ins. Default=B (FUTURE)")

        # ===Parent parser for the two variant set building parsers=================================
        validation_parent = subparsers.add_parser(name="_parentvalset",
                                                  formatter_class=CustomHelp,
                                                  add_help=False)
        # Acceptor FastQ args.
        fq1_arg = validation_parent.add_mutually_exclusive_group(required=True)
        fq1_arg.add_argument("-1", "--acceptor-fq-r1", nargs="+", dest="acceptor_fq_1s",
                             type=self.is_existing_file, metavar=("<fastqR1>", "<fastqR1_2>"),
                             help="Acceptor FastQ R1 file(s)")
        fq1_arg.add_argument("-1L", "--acceptor-fq-r1-list", dest="acceptor_fq_1s",
                             type=self.are_existing_files, metavar="<file>",
                             help="Acceptor FastQ R1 files listed per line in <file>.")
        fq2_arg = validation_parent.add_mutually_exclusive_group(required=True)
        fq2_arg.add_argument("-2", "--acceptor-fq-r2", nargs="+", dest="acceptor_fq_2s",
                             type=self.is_existing_file, metavar=("<fastqR2>", "<fastqR2_2>"),
                             help="Acceptor FastQ R2 file(s)")
        fq2_arg.add_argument("-2L", "--acceptor-fq-r2-list", dest="acceptor_fq_2s",
                             type=self.are_existing_files, metavar="<file>",
                             help="Acceptor FastQ R2 files listed per line in <file>.")
        # Optionals.
        validation_parent.add_argument("--fastq-out", metavar="<prefix>",
                                       default="VaSe_" + str(datetime.date.today()),
                                       help=("Prefix name for output FastQ files. Lane and "
                                             "pair number are appended automatically."))
        validation_parent.add_argument("--seed", default=2,
                                       type=int, metavar="<int>",
                                       help="Random seed used to randomly distribute spike-in reads. (Default='VaSe_<date>'")
        validation_parent.add_argument("-av", "--acceptor-vcf",
                                       type=self.is_variant_file, metavar="<vcf>",
                                       help="Acceptor VCF file, used to make hybrid validation VCF.")

        # ===Equivalent to AC and AB modes==========================================================
        parser_assemble = subparsers.add_parser(name="AssembleValidationSet",
                                                formatter_class=CustomHelp,
                                                parents=[univ_parent, validation_parent],
                                                help="Make validation set using acceptor FastQs and outputs from BuildSpikeIns.")
        # Varcons xor varcon list.
        vacon_arg = parser_assemble.add_mutually_exclusive_group(required=True)
        vacon_arg.add_argument("-c", "--varcon", dest="varcons_in",
                               type=self.is_existing_file, metavar="<varconfile>",
                               help="Pre-made variant context file.")
        # vacon_arg.add_argument("-c", "--varcon", nargs="+", dest="varcons_in",
        #                        type=self.is_existing_file, metavar=("<varconfile>", "<varconfile2>"),
        #                        help="Pre-made variant context file(s).")
        # vacon_arg.add_argument("-cL", "--varcon-list", dest="varcons_in",
        #                        type=self.are_existing_files, metavar="<file>",
        #                        help="Pre-made variant context files listed per line in <file>.")
        # Spike-in read files.
        spike_read_args = parser_assemble.add_mutually_exclusive_group(required=True)
        spike_read_args.add_argument("-kb", "--spike-in-bam", nargs="+", dest="spike_in_bams",
                                     type=self.is_alignment_file, metavar=("<bam>", "<bam2>"),
                                     help="Pre-built spike-in BAM file(s).")
        spike_read_args.add_argument("-kbL", "--spike-in-bam-list", dest="spike_in_bams",
                                     type=self.are_alignment_files, metavar="<file>",
                                     help="Pre-built spike-in BAM files listed per line in <file>.")
        spike_read_args.add_argument("-kfq", "--spike-in-fastq-list", dest="spike_in_fastqs",
                                     type=self.are_existing_file_pairs, metavar="<file>",
                                     help=("Pre-built spike-in FastQ files with pairs "
                                           "listed tab-separated per line in <file>."))
        # Spike-in VCF files (optional).
        spike_vcf_args = parser_assemble.add_mutually_exclusive_group()
        spike_vcf_args.add_argument("-kv", "--spike-in-vcf", nargs="+", dest="spike_in_vcfs",
                                    type=self.is_variant_file, metavar=("<vcf>", "<vcf2>"),
                                    help="Pre-built spike-in VCF file(s).")
        spike_vcf_args.add_argument("-kvL", "--spike-in-vcf-list", dest="spike_in_vcfs",
                                    type=self.are_variant_files, metavar="<file>",
                                    help="Pre-built spike-in VCF files listed per line in <file>.")

        # ===Equivalent to F mode===================================================================
        parser_full = subparsers.add_parser(name="BuildValidationSet",
                                            formatter_class=CustomHelp,
                                            parents=[univ_parent, context_parent, validation_parent],
                                            help=("Make validation set, start to finish. "
                                                  "Equivalent to BuildSpikeIns + AssembleValidationSet, "
                                                  "without intermediary, reusable spike-in files."))
        parser_full.add_argument("-a", "--acceptor-bam", required=True,
                                 type=self.is_alignment_file, metavar="<bam>",
                                 help="Acceptor BAM or CRAM file.")

    @classmethod
    def is_alignment_file(cls, file):
        """Check if path points to BAM file."""
        cls.is_existing_file(file)
        type_check = subprocess.run(
            ["file", "-b", "-z", file],
            stdout=subprocess.PIPE, check=True
            ).stdout.decode()
        if not ("BAM" in type_check or "CRAM" in type_check):
            raise argparse.ArgumentTypeError(f"File {file} is not a supported"
                                             " alignment file type.")
        with pysam.AlignmentFile(file) as infile:
            if not infile.check_index():
                raise argparse.ArgumentTypeError(f"No index found for {file}")
        return file

    @classmethod
    def are_alignment_files(cls, listfile):
        """Read list of paths and check if each is a BAM file."""
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_alignment_file(test_file)
        return files

    @classmethod
    def is_variant_file(cls, file):
        """Check if path points to VCF file."""
        cls.is_existing_file(file)
        type_check = subprocess.run(
            ["file", "-b", "-z", file],
            stdout=subprocess.PIPE, check=True
            ).stdout.decode()
        if not ("VCF" in type_check or "BCF" in type_check):
            raise argparse.ArgumentTypeError(f"File {file} is not a supported"
                                             " variant call file type.")
        with pysam.VariantFile(file) as infile:
            if not infile.index:
                raise argparse.ArgumentTypeError(f"No index found for {file}")
        return file

    @classmethod
    def are_variant_files(cls, listfile):
        """Read list of paths and check if each is a VCF file."""
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_variant_file(test_file)
        return files

    @staticmethod
    def is_existing_file(file):
        """Check if path points to existing file."""
        if not os.path.isfile(file):
            raise argparse.ArgumentTypeError(f"File {file} does not exist.")
        return file

    @staticmethod
    def is_not_existing_file(file):
        """Check if path does not point to existing file."""
        if os.path.isfile(file):
            raise argparse.ArgumentTypeError(f"File {file} already exists.")
        return file

    @classmethod
    def are_existing_files(cls, listfile):
        """Read list of paths and check if each is an existing file."""
        cls.is_existing_file(listfile)
        with open(listfile) as infile:
            files = infile.readlines()
        files = [line.strip() for line in files if not line.startswith("#")]
        for test_file in files:
            cls.is_existing_file(test_file)
        return files

    @classmethod
    def are_existing_file_pairs(cls, fastq_pair_file):
        """Read list of path pairs and check if each is an existing file."""
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
        """Check if dir exists and has write permission."""
        realpath = os.path.abspath(directory)
        if not (os.path.isdir(realpath) and os.access(realpath, os.X_OK | os.W_OK)):
            raise argparse.ArgumentTypeError(f"Directory {directory} does not exist "
                                             "or is not writeable.")
        return realpath + "/"

    @classmethod
    def is_valid_filter_file(cls, filter_file):
        """Check if file has required header."""
        cls.is_existing_file(filter_file)
        reqd_headers = ["sample", "chrom", "pos", "ref", "alt"]
        with open(filter_file) as infile:
            headers = next(infile)
        headers = [x.lower() for x in headers.strip().split("\t")]
        if headers[:5] != reqd_headers:
            raise argparse.ArgumentTypeError(f"Filter file {filter_file} header is not "
                                             "formatted properly. See help.")
        return filter_file

    @staticmethod
    def is_valid_filter_format(filter_arg):
        """Check if filter argument is formatted correctly."""
        formatted = filter_arg.split(":", 1)
        try:
            formatted[1] = [val for val in formatted[1].split(",") if val]
            cleanup = []
            for val in formatted[1]:
                if val not in cleanup:
                    cleanup.append(val)
            formatted[1] = cleanup
            formatted[0] = formatted[0].lower()
        except IndexError:
            raise argparse.ArgumentTypeError(f"Incorrect filter format: '{filter_arg}'")
        return formatted
