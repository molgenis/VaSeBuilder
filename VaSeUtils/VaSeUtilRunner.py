import os
from VariantContextFile import VariantContextFile
from VaSeUtilHelper import VaSeUtilHelper


class VaSeUtilRunner:
    def __init__(self):
        self.vuh = VaSeUtilHelper()

    def check_donor_files_exist(self, donorlistfile, samplefilter=None):
        """Checks whether the donor files used in a VaSeBuilder run still exist.

        Parameters
        ----------
        donorlistfile : str
            File with list of used donor files n VaSeBuilder run
        samplefilter : list of str
            Filter with sample names
        """
        missing_count = 0
        total_count = 0
        try:
            with open(donorlistfile, "r") as dlfile:
                next(dlfile)    # Skip the header line of the file
                for fileline in dlfile:
                    filelinedata = fileline.strip().split("\t")
                    donorfiles = filelinedata[1].split(",")

                    # Check whether the used donor files exist
                    if self.vuh.passes_filter(filelinedata[0], samplefilter):
                        for dfile in donorfiles:
                            total_count += 1
                            if not os.path.isfile(dfile):
                                missing_count += 1
                                print(f"Donor file {dfile} for sample {filelinedata[0]} could not be found. Maybe it "
                                      "has been moved, renamed or deleted?")
            print(f"Found {total_count - missing_count}/{total_count} donor files.")
        except IOError:
            print(f"Could not open {donorlistfile}")
        finally:
            print("Finished running VaSe util CheckDonorFilesExist")

    def combine_variant_context_file(self):
        print("aap")

    def generate_config_file_from_command(self, vasebuilder_command, outputpath):
        """Constructs and writes a VaSeBuilder config file based on a provided

        Parameters
        ----------
        vasebuilder_command : str
            VaSeBuilder command to construct
        outputpath : str
            Path and name to write the constructed config file
        """
        flag_split = [x for x in vasebuilder_command.split("-") if x != ""]
        for flag_value in flag_split:


    def log_info(self, vaselogloc, logfilter):
        """Displays log entries satisfying the set log level filter.

        Parameters
        ----------
        vaselogloc:
            Path to VaSeBuilder log file
        logfilter: list of str
            Filter with log level(s)
        """
        try:
            with open(vaselogloc, 'r') as vaselogfile:
                for fileline in vaselogfile:
                    fileline_elements = fileline.split("\t")
                    if len(fileline_elements) >= 3:
                        if self.vuh.passes_filter(fileline_elements[2], logfilter):
                            print(fileline.strip())
        except IOError as ioe:
            print(f"Could not open log file")

    def subset_acceptor_vcf(self):
        print("aap")

    def subset_vcf_by_variant_contexts(self, variantcontextfile, vcffile_list):
        """Subsets a filter of VCF files using a variant context files.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Variant context file to use as filter
        vcffile_list : list of str
            List of VCF files to filter
        """
        varconfile = VariantContextFile(variantcontextfile)
        for vcffile in vcffile_list:
            self.variantcontext_vcf_subsetting(varconfile, vcffile)

    def subset_vcf_by_variant_list(self, variantlistloc, vcffile_list):
        """Subsets one or more VCF files

        Parameters
        ----------
        variantlistloc : str
            Path to the variant list file to use as filter
        vcffile_list : list of str
            Paths to VCF files to subset
        """
        variantlist_filter = self.vuh.read_variant_list(variantlistloc)

        # Iterate over the VCF files to subset with the variant list
        for vcffile in vcffile_list:
            self.variantlist_vcf_subsetting(vcffile, variantlist_filter)

    def variantcontext_vcf_subsetting(self, varconfile, vcffileloc):
        """Subsets a single VCF file. The filtered entries are written to a new file.

        varconfile : VariantContextFile
            Variant contexts to use as filter
        vcffileloc : str
            Path to the VCF file to filter
        """
        tmp_data = vcffileloc.split(".")
        out_name = f"{tmp_data[0]}.varconfiltered." + ".".join(tmp_data[1:])
        try:
            vcf_file = pysam.VariantFile(vcffileloc, "r")
            filtered_vcf_file = pysam.VariantFile(out_name, "w", header=vcf_file.header)

            # Iterate over the VCF file that to be filtered nd write variants overlapping with a context
            for vcfvariant in vcf_file.fetch():
                varianttype = self.vuh.determine_variant_type(vcfvariant.ref, vcfvariant.alts)
                # max_ref = max(len(x) for x in vcfvariant.ref.split(","))
                # max_alt = max(len(x) for x in vcfvariant.alts)
                # var_endpos = max([max_ref, max_alt])
                if varconfile.variant_is_in_context(varianttype, vcfvariant.chrom, vcfvariant.start-1,
                                                    vcfvariant.stop+1):
                    filtered_vcf_file.write(vcfvariant)

            filtered_vcf_file.close()
            vcf_file.close()
        except IOError:
            print(f"Could not process variant file {vcffileloc}")

    def variantlist_vcf_subsetting(self, vcffileloc, variantlist):
        """Subsets a specified VCF file using a read variant list.

        Parameters
        ----------
        vcffileloc : str
            Path to VCF file to subset
        variantlist : dict
            Variant data to filter VCF file with
        """
        tmp_data = vcffileloc.split(".")
        out_name = f"{tmp_data[0]}.variantlist_filtered." + ".".join(tmp_data[1:])
        try:
            vcf_file = pysam.VariantFile(vcffileloc, "r")
            filtered_vcf_file = pysam.VariantFile(out_name, "w", header=vcf_file.header)

            # Iterate over the variants in the VCF file
            for vcfvariant in vcf_file.fetch():
                vcfvariant_id = f"{vcfvariant.chrom}_{vcfvariant.pos}"
                if vcfvariant_id in variantlist:
                    if self.vcf_variant_in_variantlist(vcfvariant.ref, vcfvariant.alts, variantlist[vcfvariant_id]):
                        filtered_vcf_file.write(vcfvariant)

            filtered_vcf_file.close()
            vcf_file.close()
        except IOError:
            print(f"Could not process VCF file {vcffileloc}")

    def vcf_variant_in_variantlist(self, vcfvariantref, vcfvariantalts, variantlist_entries):
        """

        Parameters
        ----------
        vcfvariantref : str
            VCF variant reference allele(s)
        vcfvariantalts : tuple of str
            VCF variant alternative allele(s)
        variantlist_entries : : list
            VCF variant(s) on the specific genomic position

        Returns
        -------
        bool
            True if variant matches entry in variant filter list, False if not
        """
        for variantlist_entry in variantlist_entries:
            if vcfvariantref == variantlist_entry[0] and vcfvariantalts == variantlist_entry[1]:
                return True
        return False
