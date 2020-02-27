#!/usr/bin/env python
"""Main VaSe running module.

This module contains the main VaSeBuilder script to parse command line options,
import necessary modules, switch modes, etc.
"""

# Import necessary modules.
import logging
import sys
import subprocess
import uuid
import time
import pysam

# Import VaSe classes.
import argparser_beta
from sample_mapper import SampleMapper
from vasebuilder import VaSeBuilder
from variant_context_file import VariantContextFile
from inclusion_filter import InclusionFilter


class VaSe:
    def __init__(self):
        """Check both the python and pysam version.

        If the python version is not at least 3.6 or higher the program will
        not run as f-strings, that are used in this program were not available
        in older versions. The program also requires pysam 0.15 or higher as
        some pysam functions used in VaseBuilder are only available since
        pysam 0.15.

        Attributes
        ----------
        """
        python_major = sys.version_info[0]
        python_minor = sys.version_info[1]
        pysam_major = int(pysam.version.__version__.split(".")[0])
        pysam_minor = int(pysam.version.__version__.split(".")[1])
        try:
            file_command = subprocess.run(["file", "-v"], check=True,
                                          stdout=subprocess.PIPE,
                                          stderr=subprocess.STDOUT).stdout
            file_version = float(file_command.decode().split("\n")[0].split("-")[1])
        except (subprocess.CalledProcessError, IndexError):
            print("Unable to detect 'file' command version. File command "
                  "likely missing or old.")
            sys.exit()
        assert (python_major >= 3 and python_minor >= 6), "Python >= 3.6 required."
        assert (pysam_major >= 0 and pysam_minor >= 15), "Pysam >= 0.15 required."
        assert file_version >= 5.37, "GNU file > 5.37 required."

        parser = argparser_beta.VctorParser()
        parser.setup()
        self.args = parser.parse_args()
        if self.args.runmode is None:  # 'required' parameter for subparsers only in Py3.7+.
            parser.parse_args(["-h"])  # This will sys.exit.
        self.vase_b = VaSeBuilder(uuid.uuid4().hex)
        self.vaselogger = self.start_logger(self.args.log, self.args.debug)

    # Runs the program.
    def main(self):
        """Run VaSeBuilder and perform all the work."""
        # Write the used command to call the program to log.
        vase_called_command = " ".join(sys.argv)
        self.vaselogger.info(f"python {vase_called_command}")

        getattr(self, self.args.runmode.lower())()

        self.vaselogger.info("VaSeBuilder run completed successfully.")
        elapsed = time.strftime("%Hh:%Mm:%Ss", time.gmtime(
            time.time() - self.vase_b.creation_time.timestamp()
            ))
        self.vaselogger.info(f"Elapsed time: {elapsed}.")

    @staticmethod
    def start_logger(logloc, debug_mode=False):
        """Start and return the logger VaSe_Logger.

        The logger writes both to stdout and the specified logfile.

        Returns
        -------
        vaselogger : Logger
            Logging utility to log VaSeBuilder activity
        """
        vaselogger = logging.getLogger("VaSe_Logger")
        vaselogger.setLevel(logging.INFO)
        vaselog_format = logging.Formatter("%(asctime)s	%(name)s	%(levelname)s	%(message)s")

        # Create the log stream to log file.
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
        self.vaselogger.info("Building sample map.")
        sample_list = SampleMapper.build_sample_maps(self.args.donor_bams,
                                                     self.args.donor_vcfs)
        variantfilter = None
        if self.args.inclusion_filter is not None:
            self.vaselogger.info("Building variant filter.")
            variantfilter = InclusionFilter.read_variant_filter_file_v2(
                self.args.inclusion_filter,
                self.args.subset_filter,
                self.args.prioritization
                )

        if not self.args.varcons_in and self.args.acceptor_bam:
            self.vaselogger.info("Building variant contexts.")
            varconfile = self.vase_b.bvcs(sample_list,
                                          self.args.acceptor_bam,
                                          self.args.out_dir,
                                          self.args.reference,
                                          self.args.varcon_out)
            if varconfile is None:
                return
        elif self.args.varcons_in and not self.args.acceptor_bam:
            # TODO: Make a way to automate multiple varcon combining here. Use VaSeUtils.MergeVarcons.py?
            varconfile = VariantContextFile(self.args.varcons_in)

        if self.args.varcon_only:
            # TODO: Modify X mode to fit here using A/D/P to change what gets output.
            pass

        elif self.args.output_mode == "A":
            # TODO: Change name of method:
            self.vase_b.run_d_mode_v2(varconfile, self.args.reference,
                                      varconfile.get_donor_alignment_files(),
                                      self.args.out_dir)

        elif self.args.output_mode == "D":
            # TODO: Make a method that's like P-mode but combines each sample.
            pass

        elif self.args.output_mode == "P":
            self.vaselogger.info("Making spike-ins per variant context.")
            self.vase_b.run_p_mode_v3(sample_list, varconfile.get_donor_alignment_files(),
                                      varconfile, self.args.out_dir)

    def assemblevalidationset(self):
        # TODO: Make a way to automate multiple varcon combining here. Use VaSeUtils.MergeVarcons.py?
        varconfile = VariantContextFile(self.args.varcons_in)
        if self.args.spike_in_bams:
            self.vase_b.run_ab_mode_v2(varconfile,
                                       self.args.acceptor_fq_1s,
                                       self.args.acceptor_fq_2s,
                                       self.args.spike_in_bams,
                                       self.args.seed,
                                       self.args.out_dir + "/" + self.args.fastq_out)
        elif self.args.spike_in_fastqs:
            self.vase_b.run_ac_mode_v2(self.args.acceptor_fq_1s,
                                       self.args.acceptor_fq_2s,
                                       self.args.spike_in_fastqs,
                                       varconfile,
                                       self.args.seed,
                                       self.args.out_dir + "/" + self.args.fastq_out)

    def buildvalidationset(self):
        sample_list = SampleMapper.build_sample_maps(self.args.donor_bams,
                                                     self.args.donor_vcfs)
        variantfilter = None
        if self.args.inclusion_filter is not None:
            variantfilter = InclusionFilter.read_variant_filter_file_v2(
                self.args.inclusion_filter,
                self.args.subset_filter,
                self.args.prioritization
                )
        varconfile = self.vase_b.bvcs(sample_list,
                                      self.args.acceptor_bam,
                                      self.args.out_dir,
                                      self.args.reference,
                                      self.args.varcon_out,
                                      variantfilter)
        self.vase_b.run_f_mode(varconfile,
                               self.args.acceptor_fq_1s,
                               self.args.acceptor_fq_2s,
                               self.args.out_dir + "/" + self.args.fastq_out,
                               self.args.seed)

    # TODO: Different 'dest' values make this hard to implement now. Try to think
    # of a way to store raw commands maybe.
# =============================================================================
#     def write_config_file(self):
#         """Write a VaSeBuilder configuration file based on the provided command line parameters.
#
#         No configuration output file will be written if one was used. The output configuration file allows
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
#             self.vaselogger.warning("Could not write config file from set command line parameters")
# =============================================================================


# Run the program.
if __name__ == "__main__":
    VASE_RUN = VaSe()
    VASE_RUN.main()
