import logging
import statistics
from VariantContext import VariantContext
from ReadIdObject import ReadIdObject


class VariantContextFile:
    """The VariantContextFile saves variant contexts

    Attributes
    -------
    vaselogger
        The logger displaying information of activities
    variant_context_file_location : str
        THe location of a read variant context file
    variant_contexts : dict
        Saves variant contexts by context identifier
    variant_context_statistics
    varcon_fields : dict
        Number representation of each variant context data field
    """

    def __init__(self, fileloc=None, samplefilter=None,
                 varconfilter=None, chromfilter=None):
        self.vaselogger = logging.getLogger("VaSe_Logger")
        self.variant_context_file_location = fileloc
        self.variant_contexts = {}
        self.variant_context_statistics = None
        self.varcon_fields = {1: "variant context id",
                              2: "sample id",
                              3: "chromosome",
                              4: "origin",
                              5: "start pos",
                              6: "end pos",
                              7: "acceptor context",
                              8: "donor context",
                              9: "number of acceptor reads",
                              10: "number of donor reads",
                              11: "acceptor/donor ratio",
                              12: "acceptor read ids",
                              13: "donor read ids"}
        if fileloc is not None:
            # Read the provided variant context file with set optional
            # filters.
            self.read_variant_context_file(fileloc, samplefilter,
                                           varconfilter, chromfilter)

    # ===METHODS TO GET DATA FROM THE VARIANT CONTEXT FILE=====================
    def get_variant_contexts(self, asdict=False):
        """Returns the variant contexts.

        By default all variant contexts are gathered in a list and returned. If the asdict parameter is set t True, the
        variant contexts are returned as a dictionary.

        Returns
        -------
        list or dict
            Variant context as list by default, as dict when parameter is set to True
        """
        if asdict:
            return self.variant_contexts
        return [varcon for varcon in self.variant_contexts.values()]

    def get_variant_contexts_by_sampleid(self):
        """Gathers variant contexts, sorts them by sample identifiers and returns them.

        Returns
        -------
        dict

        """
        varcons = self.get_variant_contexts()
        return {x.get_variant_context_sample(): [y for y in varcons if y.get_variant_context_sample() ==
                                                 x.get_variant_context_sample()] for x in varcons}

    def get_number_of_contexts(self):
        """Determines and returns the number of variant contexts.

        Returns
        -------
        int
            Number of variant contexts
        """
        return len(self.variant_contexts)

    def get_variant_context_ids(self):
        """Returns the variant context identifiers.

        Returns
        -------
        list of str
            Context identifiers of all variant contexts
        """
        return [x for x in self.variant_contexts.keys()]

    def get_variant_context(self, contextid):
        """Checks whether a variant context exists and returns it if so.

        Parameters
        ----------
        contextid : str
            Identifier of context to obtain

        Returns
        -------
        VariantContext or None
            VariantContext if it exists, None if not
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid]
        return None

    def has_variant_context(self, contextid):
        """Checks and returns whether a specified variant context has a variant context.

        Parameters
        ----------
        contextid : str
            Context identifier to check

        Returns
        -------
        bool
            True if the variant context is present, False if not
        """
        return contextid in self.variant_contexts

    def get_acceptor_context(self, contextid):
        """Fetches and returns a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Context identifier of acceptor context to return

        Returns
        -------
        OverlapContext or None
            The acceptor context if it exists, None if not
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_context()
        return None

    def get_donor_context(self, contextid):
        """Fetches and returns a specified donor context.

        Parameters
        ----------
        contextid : str
            Context identifier of donor context to return

        Returns
        -------
        OverlapContext or None
            The donor context if it exists, None if not
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_context()
        return None

    def get_all_variant_context_acceptor_reads(self):
        """Gathers and returns the acceptor reads of all variant contexts

        Returns
        -------
        list of DonorBamRead
            Acceptor reads of all variant contexts
        """
        acceptorreads = []
        for varcon in self.variant_contexts.values():
            acceptorreads.extend(varcon.get_acceptor_reads())
        return acceptorreads

    # Returns all variant context donor reads.
# =============================================================================
#     def get_all_variant_context_donor_reads(self):
#         donorreads = []
#         checklist = []
#         uniqdonorreads = []
#         for varcon in self.variant_contexts.values():
#             donorreads.extend(varcon.get_donor_reads())
#         for dbr in donorreads:
#             id_pair = (dbr.get_bam_read_id(), dbr.get_bam_read_pair_number())
#             if id_pair not in checklist:
#                 uniqdonorreads.append(dbr)
#                 checklist.append(id_pair)
#         return uniqdonorreads
# =============================================================================
    def get_all_variant_context_donor_reads(self):
        """Collects and returns all variant context donor reads.

        Returns
        -------
        list of DonorBamRead
            Variant context donor reads
        """
        dbrs = []
        donorreads = []
        for varcon in self.variant_contexts.values():
            dbrs.extend(varcon.get_donor_reads())
        for dbr in dbrs:
            donorreads.append((dbr.get_bam_read_id(),
                              dbr.get_bam_read_pair_number(),
                              dbr.get_bam_read_sequence(),
                              dbr.get_bam_read_qual()))
        return list(set(donorreads))

    def get_all_variant_context_acceptor_read_ids(self):
        """Returns all variant context acceptor read identifiers.

        Returns
        -------
        acceptorreadids : list of str
            Variant context acceptor read identifiers
        """
        acceptorreadids = []
        for varcon in self.variant_contexts.values():
            acceptorreadids.extend(varcon.get_acceptor_read_ids())
        return acceptorreadids

    def get_all_variant_context_donor_read_ids(self):
        """Returns all variant context donor read identifiers.

        Returns
        -------
        donorreadids : list of str
            Variant context donor read identifiers
        """
        donorreadids = []
        for varcon in self.variant_contexts.values():
            donorreadids.extend(varcon.get_donor_read_ids())
        return donorreadids

    def get_variant_context_fields(self):
        """Returns the number representations of the variant context data.

        Returns
        -------
        self.varcon_fields : dict
            Number representation of each variant context field
        """
        return self.varcon_fields

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA=============================
    def get_variant_context_areads(self, contextid):
        """Returns the variant context acceptor reads for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of DonorBamRead
            Variant context acceptor reads, empty list of context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_reads()
        return []

    def get_variant_context_dreads(self, contextid):
        """Returns the variant context donor reads for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of DonorBamRead
            Variant context donor reads, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_reads()
        return []

    def get_acceptor_context_reads(self, contextid):
        """Returns the reads of a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Acceptor context identifier

        Returns
        -------
        list of DonorBamRead
            Acceptor context reads
        """
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_acceptor_context():
                return self.variant_contexts[contextid].get_acceptor_context_reads()
        return []

    def get_donor_context_reads(self, contextid):
        """Returns the reads of a specified donor context.

        Parameters
        ----------
        contextid : str
            Donor context identifier

        Returns
        -------
        list of DonorBamRead
            Donor context reads, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_donor_context():
                return self.variant_contexts[contextid].get_donor_context_reads()
        return []

    # ===BASIC VARIANTCONTEXTFILE METHODS======================================
    def read_variant_context_file(self, fileloc, samplefilter=None,
                                  idfilter=None, chromfilter=None):
        """Reads a provided variant context file and saves the data.

        Filter can be set for reading the variant context file. The sample filter can nbe used to specify which samples
        to save. Samples not in the samplefilter will be skipped. Similarly filter for variant contexts and chromosome
        names can be set.

        Parameters
        ----------
        fileloc : str
            Path to variant context file to read
        samplefilter : list of str
            Sample names/identifiers to include
        idfilter : list of str
            Variant contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(fileloc, "r") as vcfile:
                varcon_records = vcfile.readlines()
        except IOError as ioe:
            self.vaselogger.critical(f"Could not read varcon file {ioe.filename}")
            exit()
        varcon_records = [x.strip().split("\t")
                          for x in varcon_records if not x.startswith("#")]

        for record in varcon_records:
            if record[0] in self.variant_contexts:
                continue

            IDpass = self.passes_filter(record[0], idfilter)
            samplepass = self.passes_filter(record[1], samplefilter)
            chrompass = self.passes_filter(record[2], chromfilter)
            if not (IDpass and samplepass and chrompass):
                continue

            acceptor_reads = [ReadIdObject(readid) for readid in record[11].split(";")]
            donor_reads = [ReadIdObject(readid) for readid in record[12].split(";")]
            new_varcon = VariantContext(*record[:6], acceptor_reads, donor_reads)
            self.variant_contexts[record[0]] = new_varcon

    def read_acceptor_context_file(self, accconfileloc, samplefilter=None, contextfilter=None, chromfilter=None):
        """Reads a provided acceptor context file and adds the acceptor contexts to already existing variant contexts.

        Parameters
        ----------
        accconfileloc : str
            Path to the acceptor contexts file
        samplefilter : list of str
            Sample names/identifiers to include
        contextfilter : list of str
            Acceptor contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(accconfileloc, "r") as accconfile:
                next(accconfile)    # Skip the header line
                for fileline in accconfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")

                    samplepass = self.passes_filter(filelinedata[1], samplefilter)
                    contextpass = self.passes_filter(filelinedata[0], contextfilter)
                    chrompass = self.passes_filter(fileline[2], chromfilter)

                    if samplepass and contextpass and chrompass:
                        if filelinedata[0] in self.variant_contexts:
                            self.variant_contexts[filelinedata[0]].add_acceptor_context(filelinedata[0], filelinedata[1]
                                                                                        , filelinedata[2],
                                                                                        int(filelinedata[3]),
                                                                                        int(filelinedata[4]),
                                                                                        int(filelinedata[5]),
                                                                                        filelinedata[7].split(";"))
        except IOError:
            self.vaselogger.warning(f"Could not read acceptor context file: {accconfileloc}")

    # Reads a donor context file
    def read_donor_context_file(self, donconfileloc, samplefilter=None, contextfilter=None, chromfilter=None):
        """Reads a provided donor contexts file and adds the donor contexts to already existing variant contexts.

        Parameters
        ----------
        donconfileloc : str
            Path to the donr contexts file
        samplefilter : list of str
            Sample names/identifiers to include
        contextfilter : list of str
            Donor contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(donconfileloc, "r") as donconfile:
                next(donconfile)    # Skip the header line
                for fileline in donconfile:
                    fileline = fileline.strip()
                    filelinedata = fileline.split("\t")

                    samplepass = self.passes_filter(filelinedata[1], samplefilter)
                    contextpass = self.passes_filter(filelinedata[0], contextfilter)
                    chrompass = self.passes_filter(filelinedata[2], chromfilter)

                    if samplepass and contextpass and chrompass:
                        if filelinedata[0] in self.variant_contexts:
                            self.variant_contexts[filelinedata[0]].add_donor_context(filelinedata[0], filelinedata[1],
                                                                                     filelinedata[2],
                                                                                     int(filelinedata[3]),
                                                                                     int(filelinedata[4]),
                                                                                     int(filelinedata[5]),
                                                                                     filelinedata[7].split(";"))
        except IOError:
            self.vaselogger.warning(f"Could not read donor context file: {donconfileloc}")

    def passes_filter(self, valtocheck, filterlist):
        """Tests and returns whether a provided value is in a provided filter.

        Filters are expected to be inclusion filters, denoting a list of values that should be used/included.

        Parameters
        ----------
        valtocheck : str
            Value to check against a filter list
        filterlist : list of str
            Values to filter with

        Returns
        -------
        bool
            True if value is in filter, False otherwise
        """
        if filterlist is not None:
            return valtocheck in filterlist
        return True

    # ===VARIANT CONTEXT SET OPERATIONS (UNION, INTERSECT, DIFFERENCE)=========
    def get_variant_contexts_union(self, other_varcon_file):
        """Returns all variant context identifiers from both variant context files.

        Parameters
        ----------
        other_varcon_file : VariantContextFile
            A VariantContextFile object with VariantContexts

        Returns
        -------
        list of str
            All variant context identifiers of both VariantContextFile objects
        """
        own_varcon_ids = self.get_variant_context_ids()
        other_varcon_ids = other_varcon_file.get_variant_context_ids()
        return list(set(own_varcon_ids) | set(other_varcon_ids))

    def get_variant_contexts_intersect(self, other_varcon_file):
        """Returns variant context identifiers present in both variant context files.

        Parameters
        ----------
        other_varcon_file : VariantContextFile
            VariantContextFile with VariantContext objects

        Returns
        -------
        list of str
            Shared variant context identifiers
        """
        own_varcon_ids = self.get_variant_context_ids()
        other_varcon_ids = other_varcon_file.get_variant_context_ids()
        return list(set(own_varcon_ids) & set(other_varcon_ids))

    def get_variant_contexts_difference(self, other_varcon_file):
        """Returns the variant context identifiers not present in the other variant context file.

        Parameters
        ----------
        other_varcon_file : VariantContextFile
            VariantContextFile with variant contexts

        Returns
        -------
        list of str
            Variant context identifiers not present in the other variant context file
        """
        own_varcon_ids = self.get_variant_context_ids()
        other_varcon_ids = other_varcon_file.get_variant_context_ids()
        return list(set(own_varcon_ids) - set(other_varcon_ids))

    def get_variant_contexts_symmetric_difference(self, other_varcon_file):
        """Determines the variant context identifiers present in either one or the other variant context file.

        THe normal difference only returns the context identifiers in A but not in B. The symmetric difference returns
        the context identifiers in A but not in B as well as context identifiers in B but not in A.

        Parameters
        ----------
        other_varcon_file : VariantContextFile

        Returns
        -------
        list of str
            Variant context read identifiers
        """
        own_varcon_ids = self.get_variant_context_ids()
        other_varcon_ids = other_varcon_file.get_variant_context_ids()
        return list(set(own_varcon_ids) ^ set(other_varcon_ids))

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA BASED ON FILTERS===============.
    def get_variant_contexts2(self, aslist=False, varconfilter=None,
                              samplefilter=None, chromfilter=None):
        """Returns the variant contexts in the VariantContextFile.

        Inclusion filters can be set to subset the returned variant contexts. FIlters include sample names/identifiers,
        variant context identifiers and chromosome names. Variant contexts can be returned as list or dictionary.

        Parameters
        ----------
        aslist : bool
            Return variant context in a list
        varconfilter : list of str
            Contexts to include
        samplefilter : list of str
            Samples to include
        chromfilter : list of str

        Returns
        -------
        list or dict of VariantContext
            Returns variant contexts as list if aslist set to True, dict otherwise
        """
        if aslist:
            return [
                    x for x in self.variant_contexts.values()
                    if (self.passes_filter(x.get_variant_context_id(),
                                           varconfilter)
                        and self.passes_filter(x.get_variant_context_sample(),
                                               samplefilter)
                        and self.passes_filter(x.get_variant_context_chrom(),
                                               chromfilter))
                    ]
        return {
                k: v for k, v in self.variant_contexts
                if (self.passes_filter(k, varconfilter)
                    and self.passes_filter(v.get_variant_context_sample(),
                                           samplefilter)
                    and self.passes_filter(v.get_variant_context_chrom(),
                                           chromfilter))
                }

    # ====METHODS TO ASSESS WHETHER A VARIANT IS IN AN EXISTING CONTEXT========
    def variant_is_in_context(self, varianttype, searchchrom,
                              searchstart, searchstop):
        """Returns whether a variant is located in an already existing context.

        If the variant type is set to 'snp' the method will check and return whether the SNP overlaps with a variant
        context. If the variant is set to 'indel' the method will check and return wheterh the area from start to end
        of the indel overlaps with a variant context.

        Parameters
        ----------
        varianttype : str
            Type of variant (snp/indel)
        searchchrom : str
            Chromosome name the variant is located on
        searchstart : int
            Leftmost genomic search window position
        searchstop : int
            Rightmost genomic search window position

        Returns
        -------
        bool or None
            True if the variant overlaps with a context, False if not
        """
        if varianttype == "snp":
            return self.snp_variant_is_in_context(searchchrom, searchstart)
        if varianttype == "indel":
            return self.indel_variant_is_in_context(searchchrom, searchstart,
                                                    searchstop)
        return None

    # Determines whether an SNP variant is located in an already existing variant context.
    def snp_variant_is_in_context(self, varchrom, vcfvarpos):
        """Checks and returns whether a SNP is located in an already existing variant context.

        For each variant context it is first checked whether the chromosome name of the context and SNP are the same. If
        so, it is then checked whether the SNP genomic position is equal or larger than the context start and smaller or
        equal than the context end.

        Parameters
        ----------
        varchrom : str
            Chromosome name the SNP is located on
        vcfvarpos : int
            Genomic position of the SNP
        Returns
        -------
        bool
            True if the SNP is in a variant context, False if not
        """
        for varcon in self.variant_contexts.values():
            if varchrom == varcon.get_variant_context_chrom():
                if (vcfvarpos >= varcon.get_variant_context_start()
                   and vcfvarpos <= varcon.get_variant_context_end()):
                    return True
        return False

    def indel_variant_is_in_context(self, indelchrom, indelleftpos, indelrightpos):
        """Checks and returns whether an indel is located in an already existing variant context.

        For each variant context it is first checked whether the chromosome name of the context and indel are the same.
        If so, it is than checked whether the indel range from indel start to end overlaps with the context.

        :Parameters
        -----------
        indelchrom : str
            The chromosome name the indel is located on
        indelleftpos : int
            The leftmost genomic position of the indel
        indelrightpos : int
            The rightmost genomic position of the indel

        Returns
        -------
        bool
            True if the indel overlaps with a variant context, False if not
        """
        for varcon in self.variant_contexts.values():
            if (indelchrom == varcon.get_variant_context_chrom()):
                if (indelleftpos <= varcon.get_variant_context_start()
                   and indelrightpos >= varcon.get_variant_context_start()):
                    return True
                if (indelleftpos <= varcon.get_variant_context_end()
                   and indelrightpos >= varcon.get_variant_context_end()):
                    return True
                if (indelleftpos >= varcon.get_variant_context_start()
                   and indelrightpos <= varcon.get_variant_context_end()):
                    return True
        return False

    def context_collision(self, context_arr):
        """Checks and returns whether a potential context overlaps with an already existing variant context.

        Parameters
        ----------
        context_arr: list
            Essential context data (chrom, variant pos, start, end)

        Returns
        -------
        bool
            True if provided context overlaps with an already existing context, False if not
        """
        if f"{context_arr[0]}_{context_arr[1]}" in self.variant_contexts:
            return True
        else:
            for varcon in self.variant_contexts.values():
                if varcon.get_variant_context_chrom() == context_arr[0]:
                    if varcon.get_variant_context_start() <= context_arr[3] \
                            and context_arr[2] <= varcon.get_variant_context_end():
                        return True
            return False

    # ===METHODS TO ADD DATA/VARIANT CONTEXTS TO THE VARIANT CONTEXT FILE======
    def set_variant_context(self, varconid, varcontext):
        """Sets a provided variant context to the provided context identifier. Will overwrite the previous value for
        the provided context identifier.

        Parameters
        ----------
        varconid : str
            Identifier of the variant context to add
        varcontext: VariantContext
            The VariantContext to set
        """
        self.variant_contexts[varconid] = varcontext

    def add_variant_context(self, varconid, sampleid,
                            varconchrom, varconorigin,
                            varconstart, varconend,
                            varcon_areads, varcon_dreads,
                            acceptor_context=None, donor_context=None):
        """Constructs a VariantContext from the provided parameters and adds it.

        Parameters
        ----------
        varconid : str
            Identifier of the variant context
        sampleid : str
            Identifier of the sample
        varconchrom : str
            Chromosome name of the variant context
        varconorigin : int
            The variant position the context is based on
        varconstart : int
            Leftmost genomic position of the variant context
        varconend : int
            Rightmost genomic position of the variant context
        varcon_areads : list of DonorBamRead
            Variant context acceptor reads
        varcon_dreads : list of DonorBamRead
            Variant contect donor reads
        acceptor_context : OverlapContext
            Acceptor context belonging to the variant context
        donor_context : OverlapContext
            Donor context belonging to the variant context
        """
        varcon_obj = VariantContext(varconid, sampleid,
                                    varconchrom, varconorigin,
                                    varconstart, varconend,
                                    varcon_areads, varcon_dreads,
                                    acceptor_context, donor_context)
        self.variant_contexts[varconid] = varcon_obj

    def add_existing_variant_context(self, varconid, varconobj):
        """Add an already created variant context to the variant context file.

        The variant context to add should be a valid VariantContext object. It will be added to the variant context
        file under the context identifier.

        Parameters
        ----------
        varconid : str
            Identifier of the context
        varconobj : VariantContext
            The created variant context
        """
        if varconobj is not None:
            self.variant_contexts[varconid] = varconobj

    def set_acceptor_context(self, varconid, acceptor_context):
        """Adds an existing acceptor context to a variant context.

        Parameters
        ----------
        varconid : str
            Context identifier
        acceptor_context : OverlapContext
            The acceptor context to add
        """
        if varconid in self.variant_contexts:
            self.variant_contexts[varconid].set_acceptor_context(acceptor_context)

    def add_acceptor_context(self, contextid, sampleid,
                             contextchrom, contextorigin,
                             contextstart, contextend,
                             acceptorreads):
        """Constructs and adds an acceptor context to a specified variant context.

        Parameters
        ----------
        contextid : str
            Context identifier
        sampleid: str
            Sample name/identifier
        contextchrom: str
            Chromosome name the context is located on
        contextorigin: int
            Variant genomic position the context is based on
        contextstart: int
            Leftmost genomic position of the context
        contextend: int
            Rightmost genomic position of the context
        acceptorreads: list of DonorBamRead
            List of reads associated with the context
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].add_acceptor_context(
                    contextid, sampleid,
                    contextchrom, contextorigin,
                    contextstart, contextend,
                    acceptorreads)

    def set_donor_context(self, varconid, donor_context):
        """Adds an already existing donor context to a specified variant context.

        Parameters
        ----------
        varconid : str
            Context identifier
        donor_context : OverlapContext
            Donor context to add
        """
        if varconid in self.variant_contexts:
            self.variant_contexts[varconid].set_donor_context(donor_context)

    def add_donor_context(self, contextid, sampleid,
                          contextchrom, contextorigin,
                          contextstart, contextend,
                          donorreads):
        """Constructs a donor context from the provided parameters and adds it to a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier
        sampleid : str
            Sample name/identifier
        contextchrom : str
            Chromosome name the context is located on
        contextorigin : int
            Variant position that constructed the context
        contextstart : int
            Donor context leftmost genomic position
        contextend : int
            Donor context rightmost genomic position
        donorreads : list of DonorBamRead
            Donor context reads
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].add_donor_context(
                    contextid, sampleid,
                    contextchrom, contextorigin,
                    contextstart, contextend,
                    donorreads)

    # ===METHODS TO ADD UNMAPPED MATE IDS TO ACCEPTOR, DONOR AND VARIANT CONTEXT=======
    def set_acceptor_context_unmapped_mate_ids(self, contextid, mateids):
        """Sets the read identifiers that have unmapped read mates fo a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Acceptor context to set unmapped read identifiers for
        mateids : lit of str
            Read identifiers with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_acceptor_context_unmapped_mates(mateids)

    def set_donor_context_unmapped_mate_ids(self, contextid, mateids):
        """Sets the read identifiers that have unmapped read mates for a specified donor context.

        Parameters
        ----------
        contextid : str
            Donor context to set unmapped read identifiers for
        mateids : list of str
            Donor context read identifiers with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_donor_context_unmapped_mates(mateids)

    def set_unmapped_acceptor_mate_ids(self, contextid, mateids):
        """Sets the acceptor read identifiers that have unmapped mates for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant contexct identifier to set unmapped acceptor read identifers for
        mateids : list of str
            Variant context acceptor read identifiers with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_acceptor_mate_ids(mateids)

    def set_unmapped_donor_mate_ids(self, contextid, mateids):
        """Sets the donor read identifiers that have unmapped mates for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier to set unmapped donor read identifiers for
        mateids: list of str
            Variant context donor read identifiers with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_donor_mate_ids(mateids)

    # ===METHODS TO GET UNMAPPED MATE IDS OF AN ACCEPTOR, DONOR AND VARIANT CONTEXT============
    def get_acceptor_context_unmapped_mate_ids(self, contextid):
        """Returns the acceptor context read identifiers that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Acceptor context identifier

        Returns
        -------
        list of str
            Acceptor context read identifiers, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_context_unmapped_mate_ids()
        return []

    def get_donor_context_unmapped_mate_ids(self, contextid):
        """Returns the donor context read identifiers that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Donor context identifier

        Returns
        -------
        list of str
            Donor contexts read identifiers, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_context_unmapped_mate_ids()
        return []

    def get_unmapped_acceptor_mate_ids(self, contextid):
        """Returns the variant context acceptor read identifiers that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of str
            Variant context acceptor read identiifers, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_acceptor_read_ids()
        return []

    def get_unmapped_donor_mate_ids(self, contextid):
        """Returns the variant context donor read identifiers that have unmapped read mates.

        Parameters
        ----------
        contextid : str

        Returns
        -------
        list of str
            Variant context donor read identifiers, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_donor_read_ids()
        return []

    # ===METHODS TO OBTAIN SOME STATISTICS ABOUT ALL THE CONTEXTS==============
    def get_average_variant_context_length(self):
        """Calculates and returns the mean variant context length.

        The mean length is calculated over all variant contexts in the variant context file. Acceptor and donor contexts
        associated with the variant contexts are not used.
        
        Returns
        -------
        float
            Average variant context length
        """
        return statistics.mean([varcon.get_variant_context_length()
                                for varcon in self.variant_contexts.values()])

    def get_median_variant_context_length(self):
        """Calculates and returns the median variant context length.

        The median length is calculated over all variant contexts in the variant context file. Acceptor and donor
        contexts associated with the variant contexts are not used.

        Returns
        -------
        float
            Median variant context length
        """
        return statistics.median([varcon.get_variant_context_length()
                                  for varcon in self.variant_contexts.values()])

    def get_average_variant_context_acceptor_reads(self):
        """Calculates and returns the average number of variant context acceptor reads.

        The mean number of acceptor reads is calculated over all variant contexts.

        Returns
        -------
        float
            Mean number of variant context acceptor reads
        """
        return statistics.mean([varcon.get_number_of_acceptor_reads() for varcon in self.variant_contexts.values()])

    def get_average_variant_context_donor_reads(self):
        """Calculates and returns the average number of variant context donor reads.

        The mean number of donor reads is calculated over all variant contexts.

        Returns
        -------
        float
            Mean number of variant context donor reads
        """
        return statistics.mean([varcon.get_number_of_donor_reads() for varcon in self.variant_contexts.values()])

    def get_median_variant_context_acceptor_reads(self):
        """Calculates and returns the median variant context acceptor reads.

        Returns
        -------
        float
            Median number of variant context acceptor reads
        """
        return statistics.median([varcon.get_number_of_acceptor_reads() for varcon in self.variant_contexts.values()])

    def get_median_variant_context_donor_reads(self):
        """Calculates and returns the median number of variant context donor reads.

        Returns
        -------
        float
            Median number of variant context donor reads
        """
        return statistics.median([varcon.get_number_of_donor_reads() for varcon in self.variant_contexts.values()])

    # ===METHODS TO WRITE VARIANT CONTEXT DATA TO A FILE=======================
    def write_variant_context_file(self, outfileloc, vbuuid, samplefilter=None,
                                   varconfilter=None, chromfilter=None):
        """Writes the data saved in the variant contexts to an output file.

        Parameters
        ----------
        outfileloc : str
            Path to write the variant context file to
        vbuuid : str
            VaSeBuilder unique identifier that build the VariantContextFile
        samplefilter : list of str
            Sample names/identifiers to include
        varconfilter : list of str
            Variant contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(outfileloc, "w") as varcon_outfile:
                varcon_outfile.write(f"#VBUUID: {vbuuid}\n")
                varcon_outfile.write("#ContextId\tDonorSample\tChrom\tOrigin\t"
                                     "Start\tEnd\tAcceptorContextLength\t"
                                     "DonorContextLength\tAcceptorReads\t"
                                     "DonorReads\tADratio\tAcceptorReadsIds\t"
                                     "DonorReadIds\n")
                for varcon in self.variant_contexts.values():
                    samplepass = self.passes_filter(
                            varcon.get_variant_context_sample(),
                            samplefilter
                            )
                    varconpass = self.passes_filter(
                            varcon.get_variant_context_id(),
                            varconfilter
                            )
                    chrompass = self.passes_filter(
                            varcon.get_variant_context_chrom(),
                            chromfilter
                            )
                    if samplepass and varconpass and chrompass:
                        varcon_outfile.write(varcon.to_string() + "\n")
        except IOError as ioe:
            self.vaselogger.warning("Could not write variant contexts to "
                                    f"{ioe.filename}")

    def write_acceptor_context_file(self, outfileloc, vbuuid, samplefilter=None,
                                    contextfilter=None, chromfilter=None):
        """Writes the acceptor contexts associated with variants contexts to an output file.

        Parameters
        ----------
        outfileloc: str
            Path and name to write acceptor contexts to
        samplefilter: list of str
            Sample names/identifiers to include
        contextfilter: list of str
            Acceptor contexts to include
        chromfilter: list of str
            Chromosome names to include
        """
        try:
            with open(outfileloc, "w") as varcon_outfile:
                varcon_outfile.write(f"VBUUID: {vbuuid}")
                varcon_outfile.write("#ContextId\tDonorSample\tChrom\tOrigin\t"
                                     "Start\tEnd\tNumOfReads\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    samplepass = self.passes_filter(
                            varcon.get_variant_context_sample(),
                            samplefilter
                            )
                    varconpass = self.passes_filter(
                            varcon.get_variant_context_id(),
                            contextfilter
                            )
                    chrompass = self.passes_filter(
                            varcon.get_variant_context_chrom(),
                            chromfilter
                            )
                    if samplepass and varconpass and chrompass:
                        varcon_outfile.write(varcon.get_acceptor_context().to_string() + "\n")
        except IOError as ioe:
            self.vaselogger.warning("Could not write acceptor contexts to "
                                    f"{ioe.filename}")

    def write_donor_context_file(self, outfileloc, vbuuid, samplefilter=None,
                                 contextfilter=None, chromfilter=None):
        """Writes the donor contexts associated with variant context to an output file.

        Parameters
        ----------
        outfileloc : str
            Path and name to write the donor contexts to
        samplefilter : list of str
            Sample names/identifiers to include
        contextfilter : list of str
            Donor contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(outfileloc, "w") as varcon_outfile:
                varcon_outfile.write(f"VBUUID: {vbuuid}")
                varcon_outfile.write("#ContextId\tDonorSample\tChrom\tOrigin\t"
                                     "Start\tEnd\tNumOfReads\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    samplepass = self.passes_filter(
                            varcon.get_variant_context_sample(),
                            samplefilter
                            )
                    varconpass = self.passes_filter(
                            varcon.get_variant_context_id(),
                            contextfilter
                            )
                    chrompass = self.passes_filter(
                            varcon.get_variant_context_chrom(),
                            chromfilter
                            )
                    if samplepass and varconpass and chrompass:
                        varcon_outfile.write(
                            varcon.get_donor_context().to_string() + "\n"
                                )
        except IOError as ioe:
            self.vaselogger.warning("Could not write donor contexts to "
                                    f"{ioe.filename}")

    def write_variant_context_stats(self, statsoutloc, vbuuid):
        """Writes basic statistics of the variant contexts to a specified output file.

        Basic statistics include means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc : str
            Path and name to write variant context statistics output file to
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"VBUUID: {vbuuid}")
                varcon_statsfile.write("#ContextId\tAvg_ALen\tAvg_DLen\t"
                                       "Med_ALen\tMed_DLen\tAvg_AQual\t"
                                       "Avg_DQual\tMed_AQual\tMed_DQual\t"
                                       "Avg_AMapQ\tAvg_DMapQ\tMed_AMapQ\t"
                                       "Med_DMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(varcon.to_statistics_string() + "\n")
        except IOError as ioe:
            self.vaselogger.critical("Could not write variant context "
                                     f"statistics to {statsoutloc}")

    def write_acceptor_context_stats(self, statsoutloc, vbuuid):
        """Writes basic statistics of the acceptor contexts to a specified output file.

        Basic statistics include means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc: str
            Path and name to write the acceptor context statistics file to
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"VBUUID: {vbuuid}")
                varcon_statsfile.write("#ContextId\tAvg_ReadLen\tMed_ReadLen\t"
                                       "Avg_ReadQual\tMed_ReadQual\t"
                                       "Avg_ReadMapQ\tMed_ReadMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(
                            varcon.get_acceptor_context().to_statistics_string()
                            + "\n"
                            )
        except IOError as ioe:
            self.vaselogger.critical("Could not write acceptor context "
                                     f"statistics to {statsoutloc}")

    # Writes some statistics about the acceptor and donor reads identified for each variant context.
    def write_donor_context_stats(self, statsoutloc, vbuuid):
        """Writes basic statistics of the donor contexts to a specified output file.

        Basic statistics include means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc : str
            Path and name to write the donor context statistics to
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"VBUUID: {vbuuid}")
                varcon_statsfile.write("#ContextId\tAvg_ReadLen\tMed_ReadLen\t"
                                       "Avg_ReadQual\tMed_ReadQual\t"
                                       "Avg_ReadMapQ\tMed_ReadMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(
                            varcon.get_donor_context().to_statistics_string()
                            + "\n"
                            )
        except IOError as ioe:
            self.vaselogger.critical("Coud not write donor context statistics "
                                     f"to {statsoutloc}")

    def write_left_right_positions(self, typetowrite, outfileloc, vbuuid):
        """Writes variant context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic positions of R1 reads and rightmost genomic positions for R2 reads are written to file.

        Parameters
        ----------
        typetowrite : str
            Type (acceptor/donor) of file to write
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"VBUUID: {vbuuid}")
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variant_contexts.values():
                    leftpositions, rightpositions = [], []
                    if typetowrite == "acceptor":
                        leftpositions = [
                                str(x)
                                for x in varcon.get_acceptor_read_left_positions()
                                ]
                        rightpositions = [
                                str(x)
                                for x in varcon.get_acceptor_read_right_positions()
                                ]
                    if typetowrite == "donor":
                        leftpositions = [
                                str(x)
                                for x in varcon.get_donor_read_left_positions()
                                ]
                        rightpositions = [
                                str(x)
                                for x in varcon.get_donor_read_right_positions()
                                ]
                    lrpof.write(str(varcon.get_variant_context_id()) + "\t"
                                + ",".join(leftpositions) + "\t"
                                + ",".join(rightpositions) + "\n")
        except IOError as ioe:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    def write_acceptor_left_right_positions(self, outfileloc, vbuuid):
        """Writes acceptor context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic position of R1 reads and rightmost genomic positions for R2 reads are written to file.

        Parameters
        ----------
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"VBUUID: {vbuuid}")
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variant_contexts.values():
                    leftpositions = [
                            str(x)
                            for x in varcon.get_acceptor_context_read_left_positions()
                            ]
                    rightpositions = [
                            str(x)
                            for x in varcon.get_acceptor_context_read_right_positions()
                            ]
                    lrpof.write(str(varcon.get_variant_context_id()) + "\t"
                                + ",".join(leftpositions) + "\t"
                                + ",".join(rightpositions) + "\n")
        except IOError as ioe:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    def write_donor_left_right_positions(self, outfileloc, vbuuid):
        """Writes donor context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic position of R1 reads and rightmost genomic positions for R2 reads are written to file.

        Parameters
        ----------
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"VBUUID: {vbuuid}")
                lrpof.write("#ContextId\tLeftPos\tRightPos\n")
                for varcon in self.variant_contexts.values():
                    leftpositions = [
                            str(x)
                            for x in varcon.get_donor_context_read_left_positions()
                            ]
                    rightpositions = [
                            str(x)
                            for x in varcon.get_donor_context_read_right_positions()
                            ]
                    lrpof.write(str(varcon.get_variant_context_id()) + "\t"
                                + ",".join(leftpositions) + "\t"
                                + ",".join(rightpositions) + "\n")
        except IOError:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    # Writes the identifiers of reads that have unmapped mates per
    # sample to a file.  Samples are all donors and the ?template?.
    def write_reads_with_unmapped_mate(self, typetowrite, umfileloc, vbuuid):
        """Writes variant context read identifiers with unmapped mates to a specified output file.

        Parameters
        ----------
        typetowrite : str
            Write variant context acceptor or donor read identifiers
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"VBUUID: {vbuuid}")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    if typetowrite == "acceptor":
                        umfile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_acceptor_mate_ids() + "\n")
                            )
                    if typetowrite == "donor":
                        umfile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_donor_mate_ids() + "\n")
                            )
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def write_acceptor_unmapped_mates(self, umfileloc, vbuuid):
        """Writes the acceptor context read identifiers that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"VBUUID: {vbuuid}")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    acccon = varcon.get_acceptor_context()
                    umfile.write(str(acccon.get_context_id()) + "\t"
                                 + str(acccon.get_sample_id()) + "\t"
                                 + ";".join(acccon.get_unmapped_read_mate_ids()) + "\n")
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def write_donor_unmapped_mates(self, umfileloc, vbuuid):
        """Writes the donor context read identifiers that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"VBUUID: {vbuuid}")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    doncon = varcon.get_donor_context()
                    umfile.write(str(doncon.get_context_id()) + "\t"
                                 + str(doncon.get_sample_id()) + "\t"
                                 + ";".join(doncon.get_unmapped_read_mate_ids())+"\n")
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def compare(self, othervarconfile, contextfilter=None):
        """Compares the current variant context file to another provided variant context file.

        Parameters
        ----------
        othervarconfile : VariantContextFile
            VariantContextFile to compare against
        contextfilter : list of str or None
            Contexts to compare (inclusive filter)

        Returns
        -------
        varcondiffs : dict
            Differences between the current and other VariantContextFile
        """
        varcondiffs = {}
        for contextid in self.variant_contexts:
            if self.passes_filter(contextid, contextfilter):
                diffs = self.variant_contexts[contextid].compare(othervarconfile.get_variant_context(contextid))
                varcondiffs[contextid] = diffs
        return varcondiffs
    
    def add_variant_context_file(self, variantcontextfile):
        """Adds another VariantContextFile to the existing one.

        Variant contexts from the second file overlapping with variant contexts from the current file will not be added.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Variant context file to add to the current
        """
        for contextid, varcon in variantcontextfile.get_variant_contexts(True).items():
            if contextid not in self.variant_contexts:
                if not self.context_collision(varcon.get_context()):
                    self.variant_contexts[contextid] = varcon
