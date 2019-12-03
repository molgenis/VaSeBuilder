#!/usr/bin/env python
import logging
import os


class ParamChecker:
    """Checks and saves the values of provided parameters by the user.

    Attributes
    ----------
    vaselogger : Logger
        Logs VaSeBuilder activity
    vcf_filelist : str
        Path to the variant list file
    bam_filelist : str
        Path to the alignment list file
    acceptorbam : str
        Path to the alignment file to use as acceptor
    fastq_in1 : list of str
        Paths to R1 fastq files to us as template
    fastq_in2 : list of str
        Paths to R2 fastq files to use as template
    outdir : str
        Path to directory to write output files to
    reference_file : str
        Path to the genome reference in fasta format
    fastq_out_location : str
        Path to write the validation fastq files to
    varcon_out_location : str
        Output name for the variant context file
    log_location : str
        Output path to write the log file to
    variantlist_location : str
    runmode : str
        The mode to run VaSeBuilder in
    varconin : str
        Path to a variant context input file
    donorfqlist : str
    required_mode_parameters : dict
        Contains the required parameters for each runmode
    """

    def __init__(self):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.vcf_filelist = ""
        self.bam_filelist = ""
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
        self.donorfqlist = ""
        self.random_seed = 2
        self.variant_filter = ""
        self.variant_priority = ()
        self.donor_aln_listfile = ""
        self.donor_var_listfile = ""
        self.bam_donor_list = ""
        self.required_mode_parameters = {"AB": ["runmode", "templatefq1", "templatefq2", "bamdonors", "varconin",
                                                "out"],
                                         "AC": ["runmode", "templatefq1", "templatefq2", "donorfastqs", "varconin",
                                                "out"],
                                         "D": ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"],
                                         "DC": ["runmode", "donorvcf", "donorbam", "out", "reference",
                                                "varconin"],
                                         "F": ["runmode", "donorvcf", "donorbam", "acceptorbam", "templatefq1",
                                               "templatefq2", "out", "reference"],
                                         "FC": ["runmode", "donorvcf", "donorbam", "templatefq1",
                                                "templatefq2", "out", "reference", "varconin"],
                                         "P": ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"],
                                         "PC": ["runmode", "donorvcf", "donorbam", "out", "reference",
                                                "varconin"],
                                         "X": ["runmode", "donorvcf", "donorbam", "acceptorbam", "out", "reference"]}
        self.optional_parameters = ["fastqout", "varcon", "variantlist", "seed", "variantfilter", "variantpriority"]

    def required_parameters_set(self, runmode, vase_arg_vals):
        """Checks and returns whether all required run mode parameters have been set.

        Parameters
        ----------
        runmode : str
            Runmode to check parameters for
        vase_arg_vals : dict
            Dictionary with parameters values

        Returns
        -------
        bool
            True if all parameters are correct, False if not
        """
        if runmode in self.required_mode_parameters:
            for reqparam in self.required_mode_parameters[runmode]:
                if vase_arg_vals[reqparam] is None:
                    self.vaselogger.critical(f"Required parameter {reqparam} is not set")
                    self.vaselogger.info("Make sure to set the required parameters: "
                                         f"{self.required_mode_parameters[runmode]}")
                    return False
            return True
        self.vaselogger.debug("Invalid run mode selected.")
        return False

    def optional_parameters_set(self, vase_arg_vals):
        """Checks and sets optional parameters if they have not been set.

        Parameters
        ----------
        vase_arg_vals:
            Dictionary with parameter values

        Returns
        -------
        vase_arg_vals : dict
            Updated dictionary with parameter values
        """
        for optparam in self.optional_parameters:
            if optparam not in vase_arg_vals:
                vase_arg_vals[optparam] = None
        return vase_arg_vals

    def check_log(self, logparam):
        """Checks the log parameter value and returns a proper output location to write the log file to.

        Parameters
        ----------
        logparam : str
            The parameter value for the log file

        Returns
        -------
        self.log_location : str
            Output path to write the log file to
        """
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

    def check_folders_exist(self, paramvals, file_exts):
        """Checks and returns whether folder containing specified file types exist.

        For each provided folder, it checks whether the folder contains at least one file with the specified file type.
        If so, the folder is added to a list.

        Parameters
        ----------
        paramvals : list of str
            Folders to check for containing specified file type
        file_exts : str
            Extension of file types to check for

        Returns
        -------
        existing_folders : list of str
            Folders containing files with the specified extension
        """
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

    def check_folder_contents(self, folder_to_check, file_exts):
        """Counts and returns the number of files with a specified extension in a specified folder.

        Parameters
        ----------
        folder_to_check:
            Folder to check for files with
        file_exts:
            File extension to check files on
        Returns
        -------
        vb_count : int
            Number of specified files in folder
        """
        vb_count = 0
        for vbfile in os.listdir(folder_to_check):
            if vbfile.endswith(file_exts):
                vb_count += 1
        self.vaselogger.debug(f"Folder {folder_to_check} contains {vb_count} "
                              f"{file_exts[0][1:4].upper()} files")
        return vb_count

    def check_file_exists(self, fileloc):
        """Checks and returns whether one or more files exist.

        Parameters
        ----------
        fileloc : str or list of str
            Path to file or list of paths to files

        Returns
        -------
        bool
            True if file(s) exists.
        """
        if type(fileloc) == str:
            fileloc = [fileloc]
        for this_file in fileloc:
            if not (os.path.isfile(this_file)):
                self.vaselogger.debug(f"File {this_file} does not exist.")
                return False
        self.vaselogger.debug("Files all exist.")
        return True

    def is_valid_output_location(self, outfilename):
        """Checks and returns whether the output location is valid.

        Parameters
        ----------
        outfilename : str
            Output file path to check

        Returns
        -------
        """
        return os.path.isdir(os.path.dirname(outfilename))

    def check_parameters(self, vase_arg_vals):
        """Checks and returns whether the set parameters are ok.

        For parameters with paths to input files, it is checked whether the file exists. For the output directory it is
        checked whether the folder exists. Optional parameters, if in the parameter value set, are set to a default
        value if none is given.

        Parameters
        ----------
        vase_arg_vals : dict
            Parameter values

        Returns
        -------
        bool
            True if set parameters are correct, False if not
        """
        # Loop over the provided parameters.
        for param in vase_arg_vals:
            self.vaselogger.debug(f"Checking parameter {param}")

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
                write_modes = ["A", "F", "D", "X", "P", "AB"]
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

                        # Check whether the priority filter and order are set as well.
                        if vase_arg_vals["variantfilter"] is not None and vase_arg_vals["variantpriority"] is not None:
                            self.variant_filter = vase_arg_vals["variantfilter"].title()
                            self.variant_priority = tuple([x.title() for x in vase_arg_vals["variantpriority"]])
                        else:
                            self.variant_filter = None
                            self.variant_priority = None
                    else:
                        self.vaselogger.warning("Variant list parameter used but supplied variant list file "
                                                f"{vase_arg_vals[param]} does not exist")

            # Checks if the provided donor fastq list file exists.
            if param == "donorfastqs":
                if vase_arg_vals[param] is not None:
                    if os.path.isfile(vase_arg_vals[param]):
                        self.donorfqlist = vase_arg_vals[param]
                    else:
                        self.vaselogger.warning("Donor fastqs parameter used but supplied donorfastq list file"
                                                f"{vase_arg_vals[param]} does not exist")

            # Checks if the provided seed value is an integer or float
            if param == "seed":
                if vase_arg_vals[param] is not None:
                    if type(vase_arg_vals[param]) is int or type(vase_arg_vals[param]) is float:
                        self.random_seed = vase_arg_vals[param]
                else:
                    self.random_seed = 2

            # Check that a valid pmode link file has been provided.
            if param == "pmodelink":
                if not os.path.isfile(vase_arg_vals[param]):
                    self.vaselogger.critical(f"P-Mode link file {vase_arg_vals[param]} does not exist")
                    return False

            # Check that a list file with donor alignment files is provided
            if param == "donoralignment":
                if vase_arg_vals[param] is not None:
                    if os.path.isfile(vase_arg_vals[param]):
                        self.donor_aln_listfile = vase_arg_vals[param]
                    else: self.donor_aln_listfile = ""
                else:
                    self.donor_aln_listfile = ""

            # Check that a list file with donor variant files is provided
            if param == "donorvariant":
                if vase_arg_vals[param] is not None:
                    if os.path.isfile(vase_arg_vals[param]):
                        self.donor_var_listfile = vase_arg_vals[param]
                    else:
                        self.donor_var_listfile = ""
                else:
                    self.donor_var_listfile = ""

            if param == "bamdonors":
                if vase_arg_vals[param] is not None:
                    if not os.path.isfile(vase_arg_vals[param]):
                        return False
                    self.bam_donor_list = vase_arg_vals[param]
                else:
                    return False
        return True

    # Only checks whether the required runmode parameters are ok using 'check_parameters()'
    def check_required_runmode_parameters(self, runmode, vaseargvals):
        """Checks and returns whether all required run modes parameters are set and correct.

        The parameter value set is first subsetted to only include all required and optional parameters. For run mode
        'X', the 'templatefq1' parameter is not required and would therefore be filtered out by this method even if it
        has been set.

        Parameters
        ----------
        runmode : str
            Selected run mode
        vaseargvals : dict
            Parameter values

        Returns
        -------
        bool
            True if all required parameters are ok, False if not
        """
        if runmode in self.required_mode_parameters:
            required_params = self.required_mode_parameters[runmode]

            # Filter command line parameters to those required for the run mode and add the optional parameters.
            filtered_vaseargvals = dict((k, v) for k, v in vaseargvals.items() if k in required_params)
            filtered_vaseargvals.update(dict((k, v) for k, v in vaseargvals.items() if k in self.optional_parameters))
            self.vaselogger.debug(f"Filtered ;parameters: {filtered_vaseargvals.items()}")
            return self.check_parameters(filtered_vaseargvals)
        return False

    def get_folder_name(self, foldername):
        """Returns the name of the folder for a given path.

        If the provided path is a folder the the provided value is returned. If the provided path is a file, the path to
        the parent folder is returned.

        Parameters
        ----------
        foldername : str
            Path to a valid folder

        Returns
        -------
        str
            Foldername parameter if foldername is a folder, the path to the parent folder if foldername is a file
        """
        if os.path.isfile(foldername) or (not os.path.isdir(foldername)):
            return os.path.dirname(foldername)
        return foldername

    # Returns the name of an output file (is used for parameters fastqout, varcon, donorbread and acceptorbread).
    def get_output_name(self, outfilename, defaultoutname):
        """Checks and returns the name of an output file.

        The output name is first checked whether it is a path. If so, only the last part (the filename) is returned. If
        no filename or path has been provided, the provided default output name is returned. Only the name, and not
        path, is returned as the output folder to write to is determined by the 'out' parameter.

        Parameters
        ----------
        outfilename : str
            The output file name
        defaultoutname : str
            The default output name to provide if non has been set/given

        Returns
        -------
        str
            Returns the output name
        """
        if outfilename is not None:
            if "/" in outfilename:
                return outfilename.split("/")[-1]
            elif "\\" in outfilename:
                return outfilename.split("\\")[-1]
            return outfilename
        return defaultoutname

    def get_valid_vcf_filelist(self):
        """Returns the location of the text file containing a list of donor variant files to use.

        Returns
        -------
        self.vcf_filelist : str
            Path to the variant list file
        """
        return self.vcf_filelist

    def get_valid_bam_filelist(self):
        """Returns the location of the text file containing a list of donor alignment files to use.

        Returns
        -------
        self.bamfile_list : str
            Path to the alignment list file
        """
        return self.bam_filelist

    def get_acceptor_bam(self):
        """Returns the location of the BAM/CRAM file that will be used as the acceptor.

        Returns
        -------
        self.acceptorbam : str
            Path to the alignment file to use as acceptor
        """
        return self.acceptorbam

    def get_first_fastq_in_location(self):
        """Returns the paths to R1 fastq files that serve as the template for the validation set.

        Returns
        -------
        self.fastq_in1 : list of str
            Paths to R1 template fastq files
        """
        return self.fastq_in1

    def get_second_fastq_in_location(self):
        """Returns the paths to R2 fastq files that serve as the template for the validation set.

        Returns
        -------
        self.fastq_in2 : list of str
            Paths to R2 template fastq files
        """
        return self.fastq_in2

    def get_fastq_in_locations(self):
        """

        Returns
        -------
        list
        """
        return [self.fastq_in1, self.fastq_in2]

    def get_out_dir_location(self):
        """Returns the directory to write the output files. Adds '/' if not present in the ouput directory location.

        Returns
        -------
        str
            Path to the output directory
        """
        if self.outdir.endswith("/"):
            return self.outdir
        return self.outdir + "/"

    def get_reference_file_location(self):
        """Returns the location of the genomic fasta reference file.

        Returns
        -------
        self.reference_file : str
            Path to the genomic reference fasta file
        """
        return self.reference_file

    # Returns the location of the FastQ file that will be produced by VaSeBuilder.
    def get_fastq_out_location(self):
        """Returns the fastq out location and suffix.

        Before returning the path and output suffix, it is checked whether the path to the outdir location ends with a
        '/'. If not this is added. The path and suffix is not a full output path. The full output path is constructed
        when writing the new validation fastq files.

        Returns
        -------
        str
            Path and suffix for the fastq out files
        """
        if self.outdir.endswith("/"):
            return self.outdir + self.fastq_out_location
        return self.outdir + "/" + self.fastq_out_location

    # Returns the location of file that will contain the variants and their context start and stops.
    def get_variant_context_out_location(self):
        """Constructs and returns the output path for the variant context file.

        Returns
        -------
        str
            Path to the variant context output location
        """
        if self.outdir.endswith("/"):
            return self.outdir + self.varcon_out_location
        return self.outdir + "/" + self.varcon_out_location

    def get_log_file_location(self):
        """Returns the location the log file will be written to.

        Returns
        -------
        self.log_location : str
            Path to write the log file to
        """
        return self.log_location

    def get_variant_list_location(self):
        """Returns the path to the file containing variants to use.

        Returns
        -------
        self.variantlist_location : str
            Path to the variant list file
        """
        return self.variantlist_location

    def get_donorfqlist(self):
        """Returns the location of the list file containing donor fastq files.

        Returns
        -------
        self.donorfqlist : str
            Path to the donor fastq list file
        """
        return self.donorfqlist

    def get_runmode(self):
        """Returns the value of the run mode parameter.

        Returns
        -------
        self.runmode : str
            The set runmode
        """
        return self.runmode

    def get_random_seed_value(self):
        """Returns the random seed value.

        Returns
        -------
        self.random_seed : int or float
            Value to use as seed for semi random distribution of donor reads of validation fastq files
        """
        return self.random_seed

    def get_variantcontext_infile(self):
        """Returns path to the variant context input file.

        Returns
        -------
        self.varconin : str
            Path to the variant context input file
        """
        return self.varconin

    def get_variant_filter(self):
        """

        Returns
        -------
        self.variant_filter : str
            Set variant priority filter
        """
        return self.variant_filter

    def get_variant_priority(self):
        """Returns the variant priority list

        Returns
        -------
        self.variant_priority
            List of variant priority labels in order
        """
        return self.variant_priority

    def get_required_runmode_parameters(self, runmode):
        """Returns the required parameters for a specified run mode.

        Parameters
        ----------
        runmode : str
            The specified run mode

        Returns
        -------
        list of str
            List with the names of the required parameters
        """
        if runmode.upper() in self.required_mode_parameters:
            return self.required_mode_parameters[runmode.upper()]
        return []

    def get_optional_parameters(self):
        """Returns the optional parameters

        Returns
        -------
        optional_parameters : list of str
            The list of optional program parameters
        """
        return self.optional_parameters

    def filter_provided_parameters(self, runmode, vase_parameters):
        """Filters and returns a parameter set with only the required and optional parameters

        Parameters
        ----------
        runmode : str
            The specified run mode
        vase_parameters : str
            Dictionary with parameter values (not used)

        Returns
        -------
        filtered_parameter_set : dict
            Parameter setf
        """
        filtered_parameter_set = {}
        runmode_reqparams = self.get_optional_parameters()
        if runmode in self.required_mode_parameters:
            runmode_reqparams.extend(self.required_mode_parameters[runmode])
        return filtered_parameter_set

    def get_donor_alignment_listfile(self):
        """Returns the path to the list file with donor alignment files.

        Returns
        -------
        self.donor_aln_files : str
            Path to list file with donor alignment files
        """
        return self.donor_aln_listfile

    def get_donor_variant_listfile(self):
        """Returns the path to the list file with donor variant files.

        Returns
        -------
        self.donor_var_files : str
            Path to list file with donor variant files
        """
        return self.donor_var_listfile

    def get_bam_donor_list(self):
        """Returns the location of the list file with donor BAM files to add to a validation set.

        Returns
        -------
        self.bam_donor_list : str
            Location of the list file with donor BAM files to add to a validation set
        """
        return self.bam_donor_list
