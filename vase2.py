#!/usr/bin/env python

# Import necessary modules.
import logging
import sys
import uuid
import argparse
import pysam
import time
from datetime import datetime

# Import VaSe classes.
from ParamChecker import ParamChecker
from VcfBamScanner2 import Sample, SampleMapper
# from VcfBamScanner import VcfBamScanner
from VaSeBuilder2 import VaSeBuilder
from VariantContextFile import VariantContextFile
from InvalidVariantFilterHeaderException import InvalidVariantFilterHeaderException
from VcfVariant import VcfVariant


class VaSe:
    def __init__(self):
        """Checks both the python and pysam version.

        If the python version is not at least 3.6 or higher the program will not run as f-strings, that are used in this
        program were not available in older versions. The program also requires pysam 0.15 or higher as some pysma
        functions used in VaseBuilder are only available since pysam 0.15.

        Attributes
        ----------
        valid_runmodes : list of str
            Valid VaSeBuilder run modes
        """
        assert (sys.version_info[0] >= 3 and sys.version_info[1] >= 6), "Please run this program in Python 3.6 or " \
                                                                        "higher"
        assert (int(pysam.version.__version__.split(".")[0]) >= 0 and int(pysam.version.__version__.split(".")[1]) >=
                15), "Please run this program with Pysam 0.15 or higher"
        self.valid_runmodes = ["AC", "D", "DC", "F", "FC", "P", "PC", "X"]

    # Runs the program.
    def main(self):
        """Runs VaSeBuilder and performs all the work.
        """
        used_config_file = False
        # Parse the command line parameters and check their validity.
        vase_arg_list = self.get_vase_parameters()
        pmc = ParamChecker()

        # Check whether a configuration file was supplied than all required command line parameters.
        if vase_arg_list["configfile"] is not None:
            used_config_file = True
            configfileloc = vase_arg_list["configfile"]
            vase_arg_list = self.read_config_file(configfileloc)    # Set the read config file as the parameter list
            # Check whether the DEBUG parameters has been set or not
            if "debug" not in vase_arg_list:
                vase_arg_list["debug"] = False
            # Check optional parameters and set those missing to None
            vase_arg_list = pmc.optional_parameters_set(vase_arg_list)

        # Start the logger and initialize this run with an ID number.
        self.vaselogger = self.start_logger(pmc, vase_arg_list["log"], vase_arg_list["debug"])

        # Write the used command to call the program to log.
        vase_called_command = " ".join(sys.argv)
        self.vaselogger.info(f"python {vase_called_command}")
        vase_b = VaSeBuilder(uuid.uuid4().hex)

        # Exit if not all of the required parameters have been set
        if not pmc.required_parameters_set(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters have been set")
            exit()

        # Exit if the supplied parameters are incorrect.
        if not pmc.check_required_runmode_parameters(vase_arg_list["runmode"], vase_arg_list):
            self.vaselogger.critical("Not all required parameters are correct. Please check log for more info.")
            exit()

        # Check if a variantfilter has been set.
        variantfilter = None
        if pmc.get_variant_list_location() != "":
            # variantfilter = self.read_variant_list(pmc.get_variant_list_location())
            variantfilter = self.read_variant_filter_file(pmc.get_variant_list_location(), pmc.get_variant_filter(),
                                                          pmc.get_variant_priority())
            # self.vaselogger.debug(f"Variant filter is {variantfilter}")

        # Check whether to write a config file from the provided command line parameters
        if not used_config_file:
            self.write_config_file(vase_arg_list)

        # Run the selected mode.
        self.run_selected_mode(pmc.get_runmode(), vase_b, pmc, variantfilter, pmc.get_random_seed_value(),
                               pmc.get_variant_filter())

        self.vaselogger.info("VaSeBuilder run completed successfully.")
        elapsed = time.strftime(
                "%Hh:%Mm:%Ss",
                time.gmtime(time.time() - vase_b.creation_time.timestamp())
                )
        self.vaselogger.info(f"Elapsed time: {elapsed}.")

    def start_logger(self, paramcheck, logloc, debug_mode=False):
        """Starts and returns the logger VaSe_Logger.

        The logger writes both to stdout and the specified logfile.

        Returns
        -------
        vaselogger : Logger
            Logging utility to log VaSeBuilder activity
        """
        vaselogger = logging.getLogger("VaSe_Logger")
        if debug_mode:
            vaselogger.setLevel(logging.DEBUG)
        else:
            vaselogger.setLevel(logging.INFO)
        vaselog_format = logging.Formatter("%(asctime)s	%(name)s	%(levelname)s	%(message)s")

        # Add the log stream to stdout.
        vase_cli_handler = logging.StreamHandler(sys.stdout)
        if debug_mode:
            vase_cli_handler.setLevel(logging.DEBUG)
        else:
            vase_cli_handler.setLevel(logging.INFO)
        vase_cli_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_cli_handler)

        # Create the log stream to log file.
        logloc = paramcheck.check_log(logloc)
        if logloc == "":
            logloc = "VaSeBuilder.log"
        vase_file_handler = logging.FileHandler(logloc)

        if debug_mode:
            vase_file_handler.setLevel(logging.DEBUG)
        else:
            vase_file_handler.setLevel(logging.INFO)
        vase_file_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_file_handler)
        return vaselogger

    def get_vase_parameters(self):
        """Creates a command line argument parser and returns the parameter values.

        Returns
        -------
        vase_args : dict
            Command line parameters and set values
        """
        # Set the VaSe parameters for the program.
        vase_argpars = argparse.ArgumentParser()
        vase_argpars.add_argument("-m", "--runmode", dest="runmode", default="F", choices=self.valid_runmodes,
                                  help="RUNMODE HELP")
        vase_argpars.add_argument("-v", "--donorvcf", dest="donorvcf",
                                  help="File containing a list of VCF/VCF.GZ/BCF files.")
        vase_argpars.add_argument("-b", "--donorbam", dest="donorbam", help="File containing a list of BAM/CRAM files.")
        vase_argpars.add_argument("-a", "--acceptorbam", dest="acceptorbam",
                                  help="BAM file for identifying acceptor reads to exclude.")
        vase_argpars.add_argument("-1", "--templatefq1", dest="templatefq1", nargs="*",
                                  help="Location and name of the first fastq in file.")
        vase_argpars.add_argument("-2", "--templatefq2", dest="templatefq2", nargs="*",
                                  help="Location and name of the second fastq in file.")
        vase_argpars.add_argument("-o", "--out", dest="out", help="Directory to write output files to.")
        vase_argpars.add_argument("-r", "--reference", dest="reference",
                                  help="Location of the reference genome. This reference genome should be used by all"
                                       "VCF/BCF and BAM/CRAM files.")
        vase_argpars.add_argument("-of", "--fastqout", dest="fastqout",
                                  help="Name for the two FastQ files to be produced.")
        vase_argpars.add_argument("-ov", "--varcon", dest="varcon",
                                  help="File name to write variants and their contexts to.")
        vase_argpars.add_argument("-l", "--log", dest="log",
                                  help="Location to write log files to (will write to working directory if not used).")
        vase_argpars.add_argument("-!", "--debug", dest="debug", action="store_true",
                                  help="Run the program in debug mode")
        vase_argpars.add_argument("-vl", "--variantlist", dest="variantlist",
                                  help="File containing a list of variants to use. Will only use these variants if "
                                       "provided. Will use all variants if no list is provided.")
        vase_argpars.add_argument("-iv", "--varconin", dest="varconin",
                                  help="Provide a Vasebuilder output variant context file to build a validation set.")
        vase_argpars.add_argument("-dq", "--donorfastqs", dest="donorfastqs",
                                  help="Location to donor fastq list file")
        vase_argpars.add_argument("-c", "--config", dest="configfile", help="Supply a config file")
        vase_argpars.add_argument("-s", "--seed", dest="seed", default=2, type=int,
                                  help="Set seed for semi randomly distributing donor reads")
        vase_argpars.add_argument("-vf", "--variantfilter", dest="variantfilter",
                                  help="")
        vase_argpars.add_argument("-vp", "--variantpriority", dest="variantpriority", nargs="*",
                                  help="Types of variant in priority.")
        vase_argpars.add_argument("-pl", "--pmodelink", dest="pmodelink",
                                  help="List file containing paths to P-mode link files")
        vase_argpars.add_argument("-da", "--donoralignment", dest="donoralignment",
                                  help="Location to list file with used donor alignment files.")
        vase_argpars.add_argument("-dv", "--donorvariant", dest="donorvariant",
                                  help="Location to list file with used donor variant files.")
        vase_argpars.add_argument("-bd", "--bamdonors", dest="bamdonors",
                                  help="List file with BAM files to add")
        vase_args = vars(vase_argpars.parse_args())
        return vase_args

    def read_variant_list(self, variantlistloc):
        """Reads a file containing genomic variants and returns them in a dictionary.

        The file containing the variant is expected to have at least three columns separated by tabs. These should be,
        in order: sample name, chromosome name, chromosomal position.

        Parameters
        ----------
        variantlistloc : str
             The location of the file containing variants

        Returns
        -------
        dict
            Read variants per sample name
        """
        variant_filter_list = {}
        try:
            with open(variantlistloc) as variantlistfile:
                next(variantlistfile)    # Skip the header line
                for fileline in variantlistfile:
                    filelinedata = fileline.strip().split("\t")
                    if filelinedata[0] not in variant_filter_list:
                        variant_filter_list[filelinedata[0]] = []
                    variant_filter_list[filelinedata[0]].append((filelinedata[1], int(filelinedata[2])))
        except IOError:
            self.vaselogger.critical(f"Could not open variant list file {variantlistloc}")
        finally:
            return variant_filter_list

    def read_config_file(self, configfileloc):
        """Reads a VaSeBuilder configuration file and returns the parameter values.

        Parameters
        ----------
        configfileloc : str
            Path to the VaSeBuilder configuration file

        Returns
        -------
        configdata : dict
            Read parameters and values
        """
        multi_value_parameters = ["templatefq1", "templatefq2", "variantpriority"]
        debug_param_vals = ["True", "1", "T"]
        configdata = {}
        try:
            with open(configfileloc, "r") as configfile:
                for fileline in configfile:
                    fileline = fileline.strip()
                    if not fileline.startswith("#"):
                        configentry = fileline.split("=")
                        if len(configentry) == 2:
                            parameter_name = configentry[0].strip().lower()
                            parameter_value = configentry[1].strip()

                            # Check whether the current parameter is a multi value parameter
                            if parameter_name in multi_value_parameters:
                                template_files = parameter_value.split(",")
                                configdata[parameter_name] = [tmplfile.strip() for tmplfile in template_files]
                            # Check if the parameter is the debug parameter
                            elif parameter_name == "debug":
                                configdata[parameter_name] = parameter_value.title() in debug_param_vals
                            # Check if the parameter is the seed parameter
                            elif parameter_name == "seed":
                                configdata[parameter_name] = int(parameter_value)
                            # What to do with the other parameters
                            else:
                                configdata[parameter_name] = parameter_value.strip()
        except IOError:
            self.vaselogger.critical(f"Could not read configuration file: {configfileloc}")
        return configdata

    # Reads a list of donor fastq files in the format (R1.fq\tR2.fq)
    def read_donor_fastq_list_file(self, donorfq_listfileloc):
        """Reads a file with a list of donor fastq files

        The donor fastq list file is expected to have two columns. The first column should contain the paths to R1 files
        and the second column paths to R2 files. For each sample two fastq files are expected,

        Parameters
        ----------
        donorfq_listfileloc : str
            Path to the donor fastq list file

        Returns
        -------
        """
        donor_fastqs = []
        try:
            with open(donorfq_listfileloc) as donorfqlistfile:
                for fileline in donorfqlistfile:
                    donor_fastqs.append(fileline.strip().split("\t"))
        except IOError:
            self.vaselogger.warning(f"Could not read donor fastq list file {donorfq_listfileloc}")
        finally:
            return donor_fastqs

    def run_selected_mode(self, runmode, vaseb, paramcheck, variantfilter, randomseed, filtercol=None):
        """Selects and runs the selected run mode.

        Depending on the specified run mode either a full validation set of fastq files or fastq files containing only
        donor data are produced. If the run mode contains a 'C' an existing variant context file will be read,
        otherwise it will be build.

        Parameters
        ----------
        runmode : str

        vaseb : VaSeBuilder
            VaSeBuilder object that will perform the actions
        paramcheck : ParamChecker
            Utility tha checks whether the parameters are ok
        variantfilter: dict
            List of variants to use as a filter
        """
        # vbscan = VcfBamScanner()
        varconfile = None    # Declare the varconfile variable so we can use it in C and no-C.

        if "X" in runmode:
            sample_list = SampleMapper.build_sample_maps(paramcheck.get_valid_bam_filelist(),
                                                         paramcheck.get_valid_vcf_filelist())
            # vcf_file_map = vbscan.scan_vcf_files(paramcheck.get_valid_vcf_filelist())
            # bam_file_map = vbscan.scan_bamcram_files(paramcheck.get_valid_bam_filelist())
            # sample_id_list = vbscan.get_complete_sample_ids()
            vaseb.run_x_mode(sample_list, paramcheck.get_acceptor_bam(),
                             paramcheck.get_reference_file_location(), paramcheck.get_out_dir_location(),
                             paramcheck.get_variant_context_out_location(), variantfilter)
        else:
            if "C" in runmode:
                # Read an existing variant context file
                varconfile = VariantContextFile(paramcheck.get_variantcontext_infile())
                if "A" in runmode:
                    donor_fastq_files = self.read_donor_fastq_list_file(paramcheck.get_donorfqlist())
                    vaseb.run_ac_mode_v2(paramcheck.get_first_fastq_in_location(),
                                         paramcheck.get_second_fastq_in_location(),
                                         donor_fastq_files, varconfile, randomseed, paramcheck.get_fastq_out_location())
                    return
                # Refetch the donor reads required when runmode (D,F,P) contains a 'C'
                # bam_file_map = vbscan.scan_bamcram_files(paramcheck.get_valid_bam_filelist())
                bam_file_map = self.read_used_donors_listfile(paramcheck.get_donor_alignment_listfile())
                vaseb.refetch_donor_reads(varconfile, bam_file_map, paramcheck.get_reference_file_location())
            elif "B" in runmode:
                varconfile = VariantContextFile(paramcheck.get_variantcontext_infile())
                if "A" in runmode:
                    # Run tbhe temporary AB mode (same aas AC but then with donor BAM files, not fastq files)
                    bam_donor_files = self.read_bam_donor_list(paramcheck.get_bam_donor_list())
                    #vaseb.run_ac_mode_v25(paramcheck.get_first_fastq_in_location(),
                    #                      paramcheck.get_second_fastq_in_location(),
                    #                      bam_donor_files, varconfile, randomseed, paramcheck.get_fastq_out_location())
                    vaseb.run_ab_mode_v2(varconfile, paramcheck.get_first_fastq_in_location(),
                                         paramcheck.get_second_fastq_in_location(), bam_donor_files,
                                         paramcheck.get_random_seed_value(), paramcheck.get_fastq_out_location())
            else:
                # Scan the variant and alignment files in the provided lists.
                sample_list = SampleMapper.build_sample_maps(paramcheck.get_valid_bam_filelist(),
                                                             paramcheck.get_valid_vcf_filelist())
# =============================================================================
#                 vcf_file_map = vbscan.scan_vcf_files(paramcheck.get_valid_vcf_filelist())
#                 bam_file_map = vbscan.scan_bamcram_files(paramcheck.get_valid_bam_filelist())
#                 sample_id_list = vbscan.get_complete_sample_ids()
# =============================================================================
                # varconfile = vaseb.build_varcon_set(sample_id_list, vcf_file_map, bam_file_map,
                #                                    paramcheck.get_acceptor_bam(), paramcheck.get_out_dir_location(),
                #                                    paramcheck.get_reference_file_location(),
                #                                    paramcheck.get_variant_context_out_location(), variantfilter)
                varconfile = vaseb.bvcs(sample_list, paramcheck.get_acceptor_bam(),
                                        paramcheck.get_out_dir_location(), paramcheck.get_reference_file_location(),
                                        paramcheck.get_variant_context_out_location(), variantfilter, filtercol)

            # Check whether contexts were created
            if varconfile is None:
                return

            # Check for modes D,F,P
            if "D" in runmode:
                vaseb.run_d_mode(varconfile, paramcheck.get_fastq_out_location())
            if "F" in runmode:
                vaseb.run_f_mode(varconfile, paramcheck.get_first_fastq_in_location(),
                                 paramcheck.get_second_fastq_in_location(), paramcheck.get_fastq_out_location(),
                                 paramcheck.get_random_seed_value())
            if "P" in runmode:
                #vaseb.run_p_mode(varconfile, paramcheck.get_out_dir_location(), paramcheck.get_fastq_out_location())
                vaseb.run_p_mode_v3(sample_list, varconfile.get_donor_alignment_files(), varconfile,
                                    paramcheck.get_out_dir_location())

    def write_config_file(self, vase_params):
        """Writes a VaSeBuilder configuration file based on the provided command line parameters.

        No configuration output file will be written if one was used. The output configuration file allows

        Parameters
        ----------
        vase_params : dict
        """
        construct_info = datetime.now()
        construct_date = construct_info.strftime("%Y%m%d")
        construct_time = construct_info.strftime("%H%M%S")
        configoutname = f"VaSe_{construct_date}_{construct_time}.cfg"
        try:
            with open(configoutname, "w") as configoutfile:
                configoutfile.write(f"#VaSe config file written on {construct_time}\n")
                for paramname, paramval in vase_params.items():
                    if vase_params[paramname] is not None:
                        if paramname == "templatefq1" or paramname == "templatefq2":
                            paramvalue = ",".join(vase_params[paramname])
                            configoutfile.write(f"{paramname.upper()}={paramvalue}\n")
                        else:
                            configoutfile.write(f"{paramname.upper()}={paramval}\n")
        except IOError:
            self.vaselogger.warning("Could not write config file from set command line parameters")

    def read_variant_filter_file(self, variant_filter_loc, priority_filter=None, priority_values=None):
        """Reads the variant filter file and saves the variants.

        The name of a column to use as a priority filter can be set. The value of this column will be saved so it can
        be used for selecting one of the overlapping variant contexts.

        Parameters
        ----------
        variant_filter_loc : str
            Path to file to use as variant filter
        priority_filter : str
            Column name to use as priority filter
        priority_values : tuple of str
            Priority values sorted in order of priority

        Returns
        -------
        variant_filter_data : dict
            Read variants from variant filter file
        """
        variant_filter_data = {}
        required_header_order = ("Sample", "Chrom", "Pos", "Ref", "Alt")

        if priority_values is not None:
            priority_values = tuple([x.title() for x in priority_values])

        try:
            with open(variant_filter_loc, "r") as variant_filter_file:
                header_line = next(variant_filter_file)
                # Check that the header is correct.
                if not self.is_valid_header(header_line, required_header_order):
                    raise InvalidVariantFilterHeaderException("Variant filter file has an incorrect header."
                                                              "Please see documentation for proper header format.")

                # Check if the set filter is in the header
                filter_column = None
                if self.filter_in_header(priority_filter, header_line):
                    self.vaselogger.debug(f"Using set priority filter {priority_filter}")
                    filter_column = self.get_filter_header_pos(priority_filter, header_line)

                # Start reading the variants from the variant filter file
                for fileline in variant_filter_file:
                    self.vaselogger.debug("Start reading variant filter file")
                    variant_data = fileline.strip().split("\t")
                    vcfvar = VcfVariant(variant_data[1], int(variant_data[2]), variant_data[3], variant_data[4])

                    # Check whether the priority filter has been set.
                    if filter_column is not None:
                        vcfvar.set_filter(priority_filter, variant_data[filter_column])
                        prlevel = self.determine_priority_index(priority_values, variant_data[filter_column])
                        vcfvar.set_priority_level(priority_filter, prlevel)

                    # Add the read variant to the filter list
                    if variant_data[0] not in variant_filter_data:
                        variant_filter_data[variant_data[0]] = []
                    variant_filter_data[variant_data[0]].append(vcfvar)
        except IOError:
            self.vaselogger.warning(f"Could not open variant filter file {variant_filter_loc}")
        except InvalidVariantFilterHeaderException:
            self.vaselogger.warning(f"Could not process variant filter file {variant_filter_loc}")
        finally:
            return variant_filter_data

    def is_valid_header(self, headerline, req_format):
        """Checks and returns whether the

        Parameters
        ----------
        headerline : str
            The header file line
        req_format : tuple of str
            Required header fields in order

        Returns
        -------
        bool
            True if header is correct, False if not
        """
        header_data = [x.title() for x in headerline.strip().split("\t")]
        for req_field in range(len(req_format)):
            if header_data[req_field].title() != req_format[req_field]:
                return False
        return True

    def filter_in_header(self, filtername, headerline):
        """Checks and returns whether

        Parameters
        ----------
        filtername : str
            Column name to use as filter
        headerline : str
            The header file line

        Returns
        -------
        bool
            True if filrter is in header, False if not
        """
        header_data = [x.title() for x in headerline.strip().split("\t")]
        if filtername is not None:
            if filtername.title() in header_data:
                return True
            return False
        return False

    def get_filter_header_pos(self, filtername, headerline):
        """Returns the index of the column to act as a priority filter.

        Parameters
        ----------
        filtername : str
            Name of the column to use as priority filter
        headerline : str
            File line containing all column headers

        Returns
        -------
        int or None
            Index of the column name if in the header, None if not
        """
        header_data = [x.title() for x in headerline.strip().split("\t")]
        if self.filter_in_header(filtername, headerline):
            return header_data.index(filtername)
        return None

    def determine_priority_index(self, priority_values, filtercolvalue):
        """Determines and returns the index of the filter column value in the ordered priority values.

        Parameters
        ----------
        priority_values : tuple of str

        filtercolvalue : str
            Value from the priority filter column

        Returns
        -------
        int or None
            The index of the value if in the priority values, None if not
        """
        if priority_values is not None and filtercolvalue is not None:
            if filtercolvalue.title() in priority_values:
                return priority_values.index(filtercolvalue.title())
        return None

    def read_used_donors_listfile(self, donor_list_file):
        """Reads a file with used alignment or variant files from a different VaSeBuilder run.

        Parameters
        ----------
        donor_list_file : str
            Path to donor list file
        """
        donor_file_map = {}
        try:
            with open(donor_list_file, "r") as dlistfile:
                next(dlistfile)
                for fileline in dlistfile:
                    filelinedata = fileline.strip().split("\t")
                    if len(fileline) == 2:
                        donor_file_map[filelinedata[0]] = filelinedata[1]
        except IOError:
            self.vaselogger.warning(f"Coud not read {donor_list_file}")
        finally:
            return donor_file_map

    def read_bam_donor_list(self, bamdlist_loc):
        """Reads the listfile with BAM files from which reads should be added to a validation set.

        Parameters
        ----------
        bamdlist_loc : str
            Path to list file with BAM files

        Returns
        -------
        dbams_to_add : list of str
        """
        dbams_to_add = []
        try:
            with open(bamdlist_loc, "r") as bamdlistfile:
                for fileline in bamdlistfile:
                    dbams_to_add.append(fileline.strip())
        except IOError:
            self.vaselogger.critical(f"")
        finally:
            return dbams_to_add


# Run the program.
if __name__ == "__main__":
    vaseb = VaSe()
    vaseb.main()
