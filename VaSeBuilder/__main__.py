#!/usr/bin/env python
"""Main VaSe running module.

This module contains the main wrapper to import necessary modules, parse
command line options, switch modes, etc.
"""

# Import necessary standard modules.
import logging
import sys
import subprocess
import uuid
import time

# Import 3rd party modules.
import pysam

# Import VaSe classes.
import argparser_beta
from sample_mapper import SampleMapper
from vasebuilder import VaSeBuilder
from variant_context_file import VariantContextFile
from inclusion_filter import InclusionFilter


class VaSe:
    """Main class to organize and run the program from the command line."""

    def __init__(self):
        """Check assertions, parse args, and initialize logger and vasebuilder.

        Checks for python >= 3.6, pysam > 0.15, and file >= 5.32. Python 3.6 is
        required for f-strings. Pysam 0.15 is required for alignment read mate
        fetching. File 5.32 is required for NGS filetype recognition.

        Attributes
        ----------
        self.args : argparse.Namespace
            Parsed arguments from custom ArgumentParser in argparser_beta module
        self.vaselogger : logging.Logger
            Logger object shared between modules
        self.vase_b : vasebuilder.VaSeBuilder
            Initialized and serialized VaSeBuilder method object
        """
        # Get python version.
        python_major = sys.version_info[0]
        python_minor = sys.version_info[1]
        # Get pysam version.
        pysam_major = int(pysam.version.__version__.split(".")[0])
        pysam_minor = int(pysam.version.__version__.split(".")[1])
        # Get system 'file' command version.
        try:
            file_command = subprocess.run(["file", "-v"], check=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT).stdout
            file_version = float(file_command.decode().split("\n")[0].split("-")[1])
        except (subprocess.CalledProcessError, IndexError):
            print("Unable to detect 'file' command version. File command "
                  "likely missing or old.")
            sys.exit()

        # Check versions and exit if below necessary versions.
        assert (python_major >= 3 and python_minor >= 6), "Python >= 3.6 required."
        assert (pysam_major >= 0 and pysam_minor >= 15), "Pysam >= 0.15 required."
        assert file_version >= 5.32, "GNU file >= 5.32 required."

        # Set up and run the argument parser.
        parser = argparser_beta.VaSeParser()
        parser.setup()
        self.args = parser.parse_args()
        if self.args.runmode is None:  # 'required' parameter for subparsers only in Py3.7+.
            print("Runmode required.")
            parser.parse_args(["-h"])  # This will sys.exit.
        if self.args.make_hash:
            try:
                import argon2
            except ModuleNotFoundError:
                print("Hashing enabled, but no argon2-cffi package found.")
        # Initialize the logger.
        self.vaselogger = self.start_logger(self.args.log, self.args.debug)
        # Initialize a VaSeBuilder instance with an ID number.
        self.vase_b = VaSeBuilder(uuid.uuid4().hex)

    def main(self):
        """Run selected VaSeBuilder methods."""
        # Log the command line used.
        vase_called_command = " ".join(sys.argv)
        self.vaselogger.info(f"python {vase_called_command}")

        # Run the selected tool.
        getattr(self, self.args.runmode.lower())()

        # Epilogue with elapsed time.
        self.vaselogger.info("VaSeBuilder run completed successfully.")
        elapsed = time.strftime("%Hh:%Mm:%Ss", time.gmtime(
            time.time() - self.vase_b.creation_time.timestamp()
            ))
        self.vaselogger.info(f"Elapsed time: {elapsed}.")

    @staticmethod
    def start_logger(logloc, debug_mode=False):
        """Start and return the logger VaSe_Logger.

        The logger writes INFO+ messages to both stdout and the specified
        logfile. DEBUG messages will print to the logfile if debug mode is
        specified.

        Returns
        -------
        vaselogger : logging.Logger
            Logging utility to log VaSeBuilder activity
        """
        # Initialize logger.
        vaselogger = logging.getLogger("VaSe_Logger")
        # Set verbosity level.
        if debug_mode:
            vaselogger.setLevel(logging.DEBUG)
        else:
            vaselogger.setLevel(logging.INFO)
        # Set log message format.
        vaselog_format = logging.Formatter("%(asctime)s	%(name)s	%(levelname)s	%(message)s")

        # Set STDOUT logging.
        vase_cli_handler = logging.StreamHandler(sys.stdout)
        vase_cli_handler.setLevel(logging.INFO)
        vase_cli_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_cli_handler)

        # Set file logging.
        if logloc is None:
            logloc = "VaSeBuilder.log"
        vase_file_handler = logging.FileHandler(logloc)
        if debug_mode:
            vase_file_handler.setLevel(logging.DEBUG)
        else:
            vase_file_handler.setLevel(logging.INFO)
        vase_file_handler.setFormatter(vaselog_format)
        vaselogger.addHandler(vase_file_handler)

        return vaselogger

    def buildspikeins(self):
        """Run BuildSpikeIns tool.

        Will produce selected outputs, such as a variant context file and/or
        spike-in BAM and VCF files, according to provided arguments.
        """
        # Check to make sure the null-run condition hasn't been set.
        if self.args.varcons_in and self.args.runmode == "V":
            self.vaselogger.warning("Requested variant contexts only, but "
                                    "also supplied pre-existing variant "
                                    "context file. No work to do. Exitting.")
            return

        # Connect BAMs and VCFs by their sample IDs.
        self.vaselogger.info("Building sample map.")
        sample_list = SampleMapper.build_sample_maps(self.args.donor_bams,
                                                     self.args.donor_vcfs,
                                                     self.args.make_hash,
                                                     self.args.hashtable)

        # Set up filter list, subsetting, and prioritization settings.
        variantfilter = None
        if self.args.inclusion_filter is not None:
            self.vaselogger.info("Building variant filter.")
            variantfilter = InclusionFilter.read_variant_filter_file_v2(
                self.args.inclusion_filter,
                self.args.subset_filter,
                self.args.prioritization
                )

        # Read pre-existing variant context file, if provided.
        if self.args.varcons_in:
            # TODO: Make a way to automate multiple varcon combining here.
            # Use VaSeUtils.MergeVarcons.py?
            varconfile = VariantContextFile(self.args.varcons_in)
            # Refetch reads and variants.
            self.vase_b.rebuild(sample_list, varconfile, self.args.reference)

        # Establish variant contexts if none provided.
        else:
            self.vaselogger.info("Building variant contexts.")
            varconfile = self.vase_b.bvcs(
                sample_list, self.args.acceptor_bam, self.args.out_dir,
                self.args.reference, self.args.varcon_out, variantfilter,
                self.args.merge
                )
            # Finish if no variant contexts were made.
            if varconfile is None:
                self.vaselogger.critical("No variant contexts built. Stopping.")
                return

        # Finish if in varcon-only mode i.e. no BAM/VCF output desired.
        if self.args.output_mode == "V":
            return

        # Write all outputs to a single BAM and single VCF file.
        if self.args.output_mode == "A":
            self.vaselogger.info("Making one combined spike-in for all contexts.")
            self.vase_b.run_a_mode_v3(sample_list, varconfile, self.args.out_dir)

        # NOT IMPLEMENTED.
        elif self.args.output_mode == "D":
            return
        #    self.vaselogger.info("Making combined spike-ins per sample.")
        #     # TODO: Make a method in between A and P that combines each SAMPLE.

        # Write each output to its own BAM and VCF file.
        elif self.args.output_mode == "P":
            self.vaselogger.info("Making spike-ins per variant context.")
            self.vase_b.run_p_mode_v3(sample_list, varconfile, self.args.out_dir)

    def assemblevalidationset(self):
        """Run AssembleValidationSet tool.

        Will produce FastQ files from spike-ins, with spike-in reads
        incorporated and acceptor reads in the corresponding loci removed.
        """
        # TODO: Make a way to automate multiple varcon combining here.
        # Use VaSeUtils.MergeVarcons.py?
        # Read existing variant context file.
        varconfile = VariantContextFile(self.args.varcons_in)

        # Write new FastQ files with donor reads added and acceptors removed.
        # Donor reads are from BAM files.
        if self.args.spike_in_bams:
            self.vase_b.run_ab_mode_v2(varconfile,
                                       self.args.acceptor_fq_1s,
                                       self.args.acceptor_fq_2s,
                                       self.args.spike_in_bams,
                                       self.args.seed,
                                       self.args.out_dir + self.args.fastq_out)
        # Donor reads are from FastQ files.
        elif self.args.spike_in_fastqs:
            self.vase_b.run_ac_mode_v2(self.args.acceptor_fq_1s,
                                       self.args.acceptor_fq_2s,
                                       self.args.spike_in_fastqs,
                                       varconfile,
                                       self.args.seed,
                                       self.args.out_dir + self.args.fastq_out)

    def buildvalidationset(self):
        """Run BuildValidationSet tool.

        Will produce a variant context file and FastQ files with spike-in
        reads incorporated and acceptor reads in the corresponding loci
        removed. Does not produced spike-in BAM and VCF files.
        """
        # Connects BAMs and VCFs by their sample IDs.
        sample_list = SampleMapper.build_sample_maps(self.args.donor_bams,
                                                     self.args.donor_vcfs,
                                                     self.args.make_hash,
                                                     self.args.hashtable)

        # Set up filter list, subsetting, and prioritization settings.
        variantfilter = None
        if self.args.inclusion_filter is not None:
            variantfilter = InclusionFilter.read_variant_filter_file_v2(
                self.args.inclusion_filter,
                self.args.subset_filter,
                self.args.prioritization
                )

        # Establish variant contexts.
        varconfile = self.vase_b.bvcs(sample_list,
                                      self.args.acceptor_bam,
                                      self.args.out_dir,
                                      self.args.reference,
                                      self.args.varcon_out,
                                      variantfilter,
                                      self.args.merge)
        # Write new FastQ files with donor reads added and acceptors removed.
        self.vase_b.run_f_mode(varconfile,
                               self.args.acceptor_fq_1s,
                               self.args.acceptor_fq_2s,
                               self.args.out_dir + self.args.fastq_out,
                               self.args.seed)

    # TODO: Different 'dest' values make this hard to implement now. Try to think
    # of a way to store raw commands maybe.
# =============================================================================
#     def write_config_file(self):
#         """Write a VaSeBuilder configuration file based on the provided command line parameters.
#
#         No configuration output file will be written if one was used. The
#         output configuration file allows
#
#         Parameters
#         ----------
#         vase_params : dict
#             Used command line parameters and values
#         """
#         construct_info = datetime.now()
#         construct_date = construct_info.strftime("%Y%m%d")
#         construct_time = construct_info.strftime("%H%M%S")
#         configoutname = f"VaSe_{construct_date}_{construct_time}.cfg"
#         try:
#             with open(configoutname, "w") as configoutfile:
#                 configoutfile.write(f"#VaSe config file written on {construct_time}\n")
#                 configoutfile.write(f"#VBUUID: {self.vase_b.creation_id}\n")
#                 for paramname, paramval in vase_params.items():
#                     if vase_params[paramname] is not None:
#                         if paramname in ["templatefq1", "templatefq2"]:
#                             paramvalue = ",".join(vase_params[paramname])
#                             configoutfile.write(f"{paramname.upper()}={paramvalue}\n")
#                         else:
#                             configoutfile.write(f"{paramname.upper()}={paramval}\n")
#         except IOError:
#             self.vaselogger.warning("Could not write config file from set "
#                                     "command line parameters")
# =============================================================================


# Run the program.
if __name__ == "__main__":
    VASE_RUN = VaSe()
    VASE_RUN.main()
