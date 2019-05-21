#!/usr/bin/env python
import logging
import os


class ParamChecker:
    # Constructor that creates two empty arrays that
    def __init__(self):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_folders = []
        self.bam_folders = []
        self.acceptorbam = ""
        self.fastq_in1 = ""
        self.fastq_in2 = ""
        self.outdir = ""
        self.fastq_out_location = ""
        self.varcon_out_location = ""
        self.log_location = ""

    # Check the logging parameter to determine where to write the
    # logfile to.
    def check_log(self, logparam):
        logloc = "VaSeBuilder.log"

        if (logparam is not None):
            # Check the location of the log file if the --log parameter
            # has been set.
            if (not(os.path.isfile(logparam)) and (logparam.endswith(".log")
                                                   or logparam.endswith(".txt"))):
                logloc = logparam

            # Check to make sure the provided --log parameter value is
            # not a directory.  (Directories could be named 'something.log').
            if (os.path.isdir(logparam)):
                logloc = logparam + "/VaSeBuilder.log"
        self.log_location = logloc
        return self.log_location

    # Checks whether provided folders exist.
    def check_folders_exist(self, paramvals, file_ext):
        existingFolders = []

        # Check whether the provided folders exist.
        for parval in paramvals:
            foldername = self.get_folder_name(parval)
            if (foldername != parval):
                self.vaselogger.warning("File instead of folder was supplied. "
                                        "Using the folder "
                                        + foldername
                                        + " as input folder")

            # Check if the supplied value is a folder or not and
            # contains any vcf/bam files.
            if (not (os.path.isdir(foldername))):
                self.vaselogger.warning("Folder "
                                        + foldername
                                        + " was not found and will therefore "
                                        "be skipped")
            else:
                if (self.check_folder_contents(foldername, file_ext) == 0):
                    self.vaselogger.warning("Folder "
                                            + foldername
                                            + " exists but contains no "
                                            + file_ext
                                            + " files")
                else:
                    self.vaselogger.info("Folder "
                                         + foldername
                                         + " will be included")
                    existingFolders.append(foldername)
        return existingFolders

    # Checks whether at least one file with a provided extension (.vcf
    # or .bam) is present.
    def check_folder_contents(self, folderToCheck, fileExt):
        vbCnt = 0
        for vbFile in os.listdir(folderToCheck):
            if (vbFile.endswith("." + fileExt)):
                vbCnt += 1
        self.vaselogger.debug("Folder "
                              + folderToCheck
                              + " contains "
                              + str(vbCnt) + " "
                              + fileExt
                              + " files")
        return vbCnt

    # Checks whether a provided file exists.
    def check_file_exists(self, fileLoc):
        if (type(fileLoc) == str):
            fileLoc = [fileLoc]
        for this_file in fileLoc:
            if (not (os.path.isfile(this_file))):
                self.vaselogger.debug("File " + this_file + " does not exist.")
                return False
        self.vaselogger.debug("Files all exist.")
        return True

    # Return the directory name of an output location.
    def is_valid_output_location(self, outFileName):
        return (os.path.isdir(os.path.dirname(outFileName)))

    # Checks whether the values of the parameters are correct (do
    # files/folders exist for example).
    # [Function should perhaps be split into smaller functions]
    def check_parameters(self, vaseArgVals):

        # Loop over the provided parameters.
        for param in vaseArgVals:

            # If the current parameter is vcfin, check whether there are
            # any valid VCF folders to use.
            if (param == 'donorvcf'):
                vcfFolders = self.check_folders_exist(vaseArgVals[param],
                                                    "vcf.gz")
                if (len(vcfFolders) == 0):
                    self.vaselogger.critical("No folders containing VCF files "
                                             "were found. Please supply "
                                             "existing folders next time :)")
                    return False
                self.vcf_folders = vcfFolders

            # If the current parameter is bamin, check whether there are
            # any valid BAM folders to use.
            if (param == 'donorbam'):
                bamFolders = self.check_folders_exist(vaseArgVals[param], "bam")
                if (len(bamFolders) == 0):
                    self.vaselogger.critical("No folders containing BAM files "
                                             "were found. Please supply "
                                             "existing folders next time :)")
                    return False
                self.bam_folders = bamFolders

            # If the current parameter is bam, check whether a valid
            # BAM file is provided.
            if (param == 'acceptorbam'):
                if (not self.check_file_exists(vaseArgVals[param])):
                    self.vaselogger.critical("No valid NIST BAM file supplied "
                                             ":(")
                    return False
                self.acceptorbam = vaseArgVals[param]

            # If the current parameter is valfastq1, check whether a
            # valid R1 fastq file is provided.
            if (param == 'templatefq1'):
                if (not self.check_file_exists(vaseArgVals[param])):
                    self.vaselogger.critical("Provided R1 FastQ input file "
                                             "does not exist")
                    return False
                self.fastq_in1 = vaseArgVals[param]

            # If the current parameter is valfastq2, check whether a
            # valid R2 fastq file is provided.
            if (param == 'templatefq2'):
                if (not self.check_file_exists(vaseArgVals[param])):
                    self.vaselogger.critical("Provided R2 FastQ input file "
                                             "does not exist")
                    return False
                self.fastq_in2 = vaseArgVals[param]

            # If the current parameter is out, check whether it is a
            # valid output location.
            if (param == 'out'):
                if (not self.is_valid_output_location(vaseArgVals[param])):
                    return False
                self.outdir = vaseArgVals[param]

            # If the current parameters is fastqout, check if a name has
            # been provided.
            if (param == 'fastqout'):
                self.fastq_out_location = self.get_output_name(vaseArgVals[param],
                                                           "VaSe")

            # If the current parameter is varcon, check whether a valid
            # output location is provided.
            if (param == 'varcon'):
                self.varcon_out_location = self.get_output_name(vaseArgVals[param],
                                                            'varcon.txt')

        # Return the lists of valid VCF and BAM folders that can be used
        # by the program.
        return True

    # Returns thename of the folder name of a parameter value (if the
    # parameter value is ).
    def get_folder_name(self, foldername):
        if (os.path.isfile(foldername) or (not os.path.isdir(foldername))):
            return os.path.dirname(foldername)
        return foldername

    # Returns the name of an output file (is used for parameters
    # fastqout, varcon, donorbread and acceptorbread).
    def get_output_name(self, outfilename, defaultoutname):
        if (outfilename is not None):
            if ("/" in outfilename):
                return outfilename.split("/")[-1]
            elif ("\\" in outfilename):
                return outfilename.split("\\")[-1]
            return outfilename
        return defaultoutname

    # Returns the list of valid VCF folders.
    def get_valid_vcf_folders(self):
        return self.vcf_folders

    # Returns the list of valid BAM folders.
    def get_valid_bam_folders(self):
        return self.bam_folders

    # Returns the location of the  NIST BAM file.
    def get_acceptor_bam(self):
        return self.acceptorbam

    # Returns the location and name of the first (R1) fastq input file.
    def get_first_fastq_in_location(self):
        return self.fastq_in1

    # Returns the location and name of the second (R2) fastq input file.
    def get_second_fastq_in_location(self):
        return self.fastq_in2

    # Returns the location(s) and names of the two (R1 and R2) fastq
    # input files.
    def get_fastq_in_locations(self):
        return [self.fastq_in1, self.fastq_in2]

    # Returns the location to write the output to.
    def get_out_dir_location(self):
        return self.outdir

    # Returns the location of the FastQ file that will be produced by
    # VaSeBuilder.
    def get_fastq_out_location(self):
        return self.outdir + "/" + self.fastq_out_location

    # Returns the location of file that will contain the variants and
    # their context start and stops.
    def get_variant_context_out_location(self):
        return self.outdir + "/" + self.varcon_out_location

    # Retuns the location to write the log file(s) to.
    def get_log_file_location(self):
        return self.log_location
