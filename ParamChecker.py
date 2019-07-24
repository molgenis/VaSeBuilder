#!/usr/bin/env python
import logging
import os


class ParamChecker:
    # Constructor that creates two empty arrays that
    def __init__(self):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_filelist = []
        self.bam_filelist = []
        self.acceptorbam = ""
        self.fastq_in1 = ""
        self.fastq_in2 = ""
        self.outdir = ""
        self.reference_file = ""
        self.fastq_out_location = ""
        self.varcon_out_location = ""
        self.log_location = ""
        self.variantlist_location = ""
        self.runmode = ""
        self.varconin = ""
        self.required_mode_parameters = {"A": ["templatefq1", "templatefq2", "donorfastqs", "varconin"],
                                         "AC": ["templatefq1", "templatefq2", "donorfastqs", "varconin"],
                                         "D": ["donorvcf", "donorbam", "acceptorbam", "out", "reference"],
                                         "DC": ["donorvcf", "donorbam", "acceptorbam", "out", "reference", "varconin"],
                                         "F": ["donorvcf", "donorbam", "acceptorbam", "templatefq1", "templatefq2",
                                               "out", "reference"],
                                         "FC": ["donorvcf", "donorbam", "acceptorbam", "templatefq1", "templatefq2",
                                                "out", "reference", "varconin"],
                                         "P": ["donorvcf", "donorbam", "acceptorbam", "out", "reference"],
                                         "PC": ["donorvcf", "donorbam", "acceptorbam", "out", "reference", "varconin"],
                                         "X": ["donorvcf", "donorbam", "acceptorbam", "out", "reference"],
                                         "XC": ["donorvcf", "donorbam", "acceptorbam", "out", "reference", "varconin"]}

    # Checks whether the required parameters for a run mode have been set.
    def required_parameters_set(self, runmode, vase_arg_vals):
        if runmode in self.required_mode_parameters:
            for reqparam in self.required_mode_parameters[runmode]:
                if vase_arg_vals[reqparam] is None:
                    self.vaselogger.debug(f"Required parameter {reqparam} is not set")
                    return False
            return True
        self.vaselogger.debug("Invalid run mode selected.")
        return False

    # Check the logging parameter to determine where to write the logfile to.
    def check_log(self, logparam):
        logloc = "VaSeBuilder.log"

        if logparam is not None:
            # Check the location of the log file if the --log parameter has been set.
            if (not(os.path.isfile(logparam)) and (logparam.endswith(".log")
                                                   or logparam.endswith(".txt"))):
                logloc = logparam

            # Check to make sure the provided --log parameter value is
            # not a directory. (Directories could be named "something.log").
            if os.path.isdir(logparam):
                logloc = logparam + "/VaSeBuilder.log"
        self.log_location = logloc
        return self.log_location

    # Checks whether provided folders exist.
    def check_folders_exist(self, paramvals, file_exts):
        existing_folders = []

        # Check whether the provided folders exist.
        for parval in paramvals:
            foldername = self.get_folder_name(parval)
            if foldername != parval:
                self.vaselogger.warning("File instead of folder was supplied. "
                                        f"Using the folder {foldername} "
                                        "as input folder")

            # Check if the supplied value is a folder or not and contains any vcf/bam files.
            if not (os.path.isdir(foldername)):
                self.vaselogger.warning(f"Folder {foldername} was not found "
                                        "and will therefore be skipped")
            else:
                if self.check_folder_contents(foldername, file_exts) == 0:
                    self.vaselogger.warning(f"Folder {foldername} exists but "
                                            f"contains no {file_exts[0][1:4].upper()} files")
                else:
                    self.vaselogger.info(f"Folder {foldername} will be "
                                         "included")
                    existing_folders.append(foldername)
        return existing_folders

    # Checks whether at least one file with a provided extension (.vcf or .bam) is present.
    def check_folder_contents(self, folder_to_check, file_exts):
        vb_count = 0
        for vbfile in os.listdir(folder_to_check):
            if vbfile.endswith(file_exts):
                vb_count += 1
        self.vaselogger.debug(f"Folder {folder_to_check} contains {vb_count} "
                              f"{file_exts[0][1:4].upper()} files")
        return vb_count

    # Checks whether a provided file exists.
    def check_file_exists(self, fileloc):
        if type(fileloc) == str:
            fileloc = [fileloc]
        for this_file in fileloc:
            if not (os.path.isfile(this_file)):
                self.vaselogger.debug(f"File {this_file} does not exist.")
                return False
        self.vaselogger.debug("Files all exist.")
        return True

    # Return the directory name of an output location.
    def is_valid_output_location(self, outfilename):
        return os.path.isdir(os.path.dirname(outfilename))

    # Checks whether the values of the parameters are correct (do files/folders exist for example).
    # [Function should perhaps be split into smaller functions]
    def check_parameters(self, vase_arg_vals):

        # Loop over the provided parameters.
        for param in vase_arg_vals:

            # If the current parameter is donorvcf, check that the file containing the list of donor VCFs exists.
            if param == "donorvcf":
                if not os.path.isfile(vase_arg_vals["donorvcf"]):
                    self.vaselogger.critical("No VCF/BCF donor list file found")
                    return False
                self.vcf_filelist = vase_arg_vals["donorvcf"]

            # If the current parameter is donorbam, check that the file containing the list of donor BAMs exists.
            if param == "donorbam":
                if not os.path.isfile(vase_arg_vals["donorbam"]):
                    self.vaselogger.critical("No BAM/CRAM donor list file found")
                    return False
                self.bam_filelist = vase_arg_vals["donorbam"]

            # If the current parameter is acceptorbam, check whether a valid BAM file is provided.
            if param == "acceptorbam":
                if not self.check_file_exists(vase_arg_vals[param]):
                    self.vaselogger.critical("No valid acceptor/template BAM file supplied :(")
                    return False
                self.acceptorbam = vase_arg_vals[param]

            # If the current parameter is valfastq1, check whether one or more valid R1 fastq files are provided.
            if param == "templatefq1":
                if not self.check_file_exists(vase_arg_vals[param]):
                    self.vaselogger.critical("No valid R1 FastQ file(s) provided")
                    return False
                self.fastq_in1 = vase_arg_vals[param]

            # If the current parameter is valfastq2, check whether one or more valid R2 fastq files are provided.
            if param == "templatefq2":
                if not self.check_file_exists(vase_arg_vals[param]):
                    self.vaselogger.critical("No valid R2 FastQ file(s) provided")
                    return False
                self.fastq_in2 = vase_arg_vals[param]

            # If the current parameter is out, check whether it is a valid output location.
            if param == "out":
                if not self.is_valid_output_location(vase_arg_vals[param]):
                    return False
                self.outdir = vase_arg_vals[param]

            # If the current parameter is reference check if the reference file exists
            if param == "reference":
                if not self.check_file_exists(vase_arg_vals[param]):
                    self.vaselogger.critical("No valid reference file supplied")
                    return False
                self.reference_file = vase_arg_vals[param]

            # If the current parameters is fastqout, check if a name has been provided.
            if param == "fastqout":
                self.fastq_out_location = self.get_output_name(vase_arg_vals[param], "VaSe")

            # If the current parameter is varcon, check whether a valid output location is provided.
            if param == "varcon":
                self.varcon_out_location = self.get_output_name(vase_arg_vals[param], "varcon.txt")

            if param == "runmode":
                write_modes = ["F", "D", "X", "P"]
                write_modes = write_modes + [x + "C" for x in write_modes]
                if vase_arg_vals[param] in write_modes:
                    self.runmode = vase_arg_vals[param]
                else:
                    self.vaselogger.critical("Invalid run-mode option. See help for more info.")
                    return False

            if param == "varconin":
                if vase_arg_vals[param] is not None:
                    if not self.check_file_exists(vase_arg_vals[param]):
                        self.vaselogger.critical("No valid variant context file supplied")
                        return False
                self.varconin = vase_arg_vals[param]

            # Checks if the provided variant list file exists
            if param == "variantlist":
                if vase_arg_vals[param] is not None:
                    if self.check_file_exists(vase_arg_vals[param]):
                        self.variantlist_location = vase_arg_vals[param]
                    else:
                        self.vaselogger.warning("Variant list parameter used but supplied variant list file "
                                                f"{vase_arg_vals[param]} does not exist")
        return True

    # Returns the name of the folder name of a parameter value (if the parameter value is ).
    def get_folder_name(self, foldername):
        if os.path.isfile(foldername) or (not os.path.isdir(foldername)):
            return os.path.dirname(foldername)
        return foldername

    # Returns the name of an output file (is used for parameters fastqout, varcon, donorbread and acceptorbread).
    def get_output_name(self, outfilename, defaultoutname):
        if outfilename is not None:
            if "/" in outfilename:
                return outfilename.split("/")[-1]
            elif "\\" in outfilename:
                return outfilename.split("\\")[-1]
            return outfilename
        return defaultoutname

    # Returns the list of valid VCF folders.
    def get_valid_vcf_filelist(self):
        return self.vcf_filelist

    # Returns the list of valid BAM folders.
    def get_valid_bam_filelist(self):
        return self.bam_filelist

    # Returns the location of the acceptor BAM file location.
    def get_acceptor_bam(self):
        return self.acceptorbam

    # Returns the location and name of the first (R1) fastq input file location.
    def get_first_fastq_in_location(self):
        return self.fastq_in1

    # Returns the location and name of the second (R2) fastq input file location.
    def get_second_fastq_in_location(self):
        return self.fastq_in2

    # Returns the location(s) and names of the two (R1 and R2) fastq input files.
    def get_fastq_in_locations(self):
        return [self.fastq_in1, self.fastq_in2]

    # Returns the location to write the output to.
    def get_out_dir_location(self):
        if self.outdir.endswith("/"):
            return self.outdir
        return self.outdir + "/"

    # Returns the location of the reference fasta file
    def get_reference_file_location(self):
        return self.reference_file

    # Returns the location of the FastQ file that will be produced by VaSeBuilder.
    def get_fastq_out_location(self):
        if self.outdir.endswith("/"):
            return self.outdir + self.fastq_out_location
        return self.outdir + "/" + self.fastq_out_location

    # Returns the location of file that will contain the variants and their context start and stops.
    def get_variant_context_out_location(self):
        if self.outdir.endswith("/"):
            return self.outdir + self.varcon_out_location
        return self.outdir + "/" + self.varcon_out_location

    # Retuns the location to write the log file(s) to.
    def get_log_file_location(self):
        return self.log_location

    # Returns the variant list location
    def get_variant_list_location(self):
        return self.variantlist_location
