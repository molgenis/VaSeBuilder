import logging

"""Contains methods that are used by multiple other classes and are therefore collected in a single helper."""


class VaSeUtilHelper:
    def __init__(self):
        self.vaseutillogger = logging.getLogger("VaSeUtil_Logger")
        self.parameter_map = {"-m": "RUNMODE",
                              "--runmode": "RUNMODE",
                              "-v": "DONORVCF",
                              "--donorvcf": "DONORVCF",
                              "-b": "DONORBAM",
                              "--donorbam": "DONORBAM",
                              "-a": "ACCEPTORBAM",
                              "--acceptorbam": "ACCEPTORBAM",
                              "-1": "TEMPLATEFQ1",
                              "--templatefq1": "TEMPLATEFQ1",
                              "-2": "TEMPLATEFQ2",
                              "--templatefq2": "TEMPLATEFQ2",
                              "-o": "OUT",
                              "--out": "OUT",
                              "-r": "REFERENCE",
                              "--reference": "REFERENCE",
                              "-of": "",
                              "--fastqout": "",
                              "-ov": "",
                              "--varcon": "",
                              "-l": "LOG",
                              "--log": "LOG",
                              "-!": "DEBUG",
                              "-vl": "VARIANTLIST",
                              "--variantlist": "VARIANTLIST",
                              "-iv": "VARCONIN",
                              "--varconin": "VARCONIN",
                              "-dq": "DONORFASTQS",
                              "--donorfastqs": "DONORFASTQS",
                              "-c": "CONFIG",
                              "--config": "CONFIG",
                              "-s": "SEED",
                              "--seed": "SEED"}

    # Returns whether something is in the filter or not
    def passes_filter(self, valtocheck, filterlist):
        if filterlist is not None:
            return valtocheck in filterlist
        return True

    # Displays BAM file information for a list of reads
    def get_bamread_info(self, readslist, searchchrom, searchstart, searchend, bamfile, readidfilter=None):
        for bamread in bamfile.fetch(searchchrom, searchstart, searchend):
            if bamread.query_name in readslist and self.passes_filter(bamread, readidfilter):
                print(bamread.to_string())

    # Reads the file with a list of used donor VCF/BAM files
    def read_donor_list_file(self, dlistfile, samplefilter=None):
        donorfiles = {}
        try:
            with open(dlistfile, 'r') as dlfile:
                next(dlfile)    # Skip the header line
                for fileline in dlfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    
                    # Check if the entry is in the set sample filter
                    if self.passes_filter(filelinedata[0], samplefilter):
                        donorfiles[filelinedata[0]] = filelinedata[1]
        except IOError as ioe:
            self.vaseutillogger.critical("Could not read donor list file")
        return donorfiles

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck').
    def read_abamreads_list_nofilter(self, acceptorreadfile):
        acceptorreads = []
        with open(acceptorreadfile, 'r') as arfile:
            next(arfile)    # Skip the header line
            for fileline in arfile:
                fileline = fileline.strip()
                filelinedata = fileline.split("\t")
                acceptorreads.extend(filelinedata[1:])
        return acceptorreads

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck').
    def read_dbamreads_list_nofilter(self, donorReadFile):
        donorReads = []
        with open(donorReadFile, 'r') as arFile:
            next(arFile)    # Skip the header line
            for fileLine in arFile:
                fileLine = fileLine.strip()
                fileLineData = fileLine.split("\t")
                donorReads.extend(fileLineData[2:])
        return donorReads

    # Reads the file containing acceptor/donor BAM reads (used for utils such as 'acceptorcheck', 'donorcheck').
    def readDBamReadsListByVarcon_noFilter(self, acceptorReadFile):
        donorReads = []
        with open(acceptorReadFile, 'r') as arFile:
            next(arFile)    # Skip the header line
            for fileLine in arFile:
                fileLine = fileLine.strip()
                fileLineData = fileLine.split("\t")
                donorReads.extend(fileLineData[2:])
        return donorReads

    # Returns the map of acceptor/donor context fields (can be used for compare results)
    def get_accdon_context_fields(self):
        accdon_fields = {1: "Context ID",
                         2: "Sample ID",
                         3: "Context chrom",
                         4: "Context origin",
                         5: "Context start",
                         6: "Context end",
                         7: "Context reads"}
        return accdon_fields

    # Returns the map of variant context fields (can be used for compare results)
    def get_variant_context_fields(self):
        varcon_fields = {}
        return varcon_fields

    # Reads a provided donor list file into a hashmap
    def read_donorlist_file(self, dllist_fileloc):
        dlist_data = {}
        try:
            with open(dllist_fileloc, "r") as dlistfile:
                next(dlistfile)    # Skip the header line
                for fileline in dlistfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")
                    dlist_data[filelinedata[0]] = filelinedata[1].split(";")
        except IOError:
            self.vaseutillogger.warning(f"Could not read the provided donor list file {dllist_fileloc}")
        return dlist_data

    # Reads the unmapped mate file
    def read_unmapped_matefile(self, umapfileloc):
        unmapped_mateids = {}
        try:
            with open(umapfileloc, 'r') as umapfile:
                next(umapfile)  # Skip the header line
                for fileline in umapfile:
                    filelinedata = fileline.strip().split("\t")
                    unmapped_mateids[filelinedata[0]] = filelinedata[2].split(";")
        except IOError:
            self.vaseutillogger.warning("Could not read acceptor unmapped mate file")

    # Returns whether the read is the first or second read in a pair.
    def get_read_pair_num(self, pysam_bamread):
        if pysam_bamread.is_read1:
            return "1"
        return "2"

    def determine_variant_type(self, vcfvariantref, vcfvariantalts):
        """Determines and returns the

        Parameters
        ----------
        vcfvariantref : str
            Variant reference allele(s)
        vcfvariantalts : tuple of str
            Variant alternative allele(s)

        Returns
        -------
        str
            Type of variant (snp/indel)
        """
        maxreflength = max([len(x) for x in vcfvariantref.split(",")])
        maxaltlength = max([len(x) for x in vcfvariantalts])

        # Check based on the reference and alternative lengths whether the variant is a SNP or indel.
        if maxreflength == 1 and maxaltlength == 1:
            return "snp"
        elif maxreflength > 1 or maxaltlength > 1:
            return "indel"
        return "?"

    def get_config_param_name(self, paramflag):
        """Returns the config parameter name for a parameter flag.

        Parmeters
        ---------
        paramflag : str
            Parameter flag to obtain config parameter

        Returns
        -------
        str
            Config parameter name if parameter flag is in map, empty string otherwise
        """
        if paramflag in self.parameter_map:
            return self.parameter_map[paramflag]
        return ""
