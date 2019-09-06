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
        Fields
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

    # Returns the variant contexts by sample identifier
    def get_variant_contexts_by_sampleid(self):
        """Gathers variant contexts, sorts them by sample identifiers and returns them.

        Returns
        -------
        dict

        """
        varcons = self.get_variant_contexts()
        return {x.get_variant_context_sample(): [y for y in varcons if y.get_variant_context_sample() ==
                                                 x.get_variant_context_sample()] for x in varcons}

    # Returns the number of contexts saved
    def get_number_of_contexts(self):
        """Determines and returns the number of variant contexts.

        Returns
        -------
        int
            Number of variant contexts
        """
        return len(self.variant_contexts)

    # Returns a list of context ids
    def get_variant_context_ids(self):
        """Returns the variant context identifiers.

        Returns
        -------
        list of str
            Context identifiers of all variant contexts
        """
        return [x for x in self.variant_contexts.keys()]

    # Returns a specified variant context.
    def get_variant_context(self, contextid):
        """Checks whether a variant context exists and returns it if so."""
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid]
        return None

    # Returns whether a variant context is present
    def has_variant_context(self, contextid):
        """Checks and returns whether a specified variant context has a variant context.

        Parameters
        ----------
        contextid : str
            Context identifier
        """
        return contextid in self.variant_contexts

    # Returns the acceptor context of the specified variant context.
    def get_acceptor_context(self, contextid):
        """Fetches and returns a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Context identifier of acceptor context to return

        Returns
        -------
        OverlapContext or None
            The acceptor context is exists, None if not
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_context()
        return None

    # Returns the donor context of the specified variant context.
    def get_donor_context(self, contextid):
        """Fetches and returns a specified donor context.

        Parameters
        ----------
        contextid : str
            Context identifier of donor context to return
        Returns
        -------
        OverlapContext or None
            The donor context if exists, None if not
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_context()
        return None

    # Returns all variant context acceptor reads.
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

    # Returns all variant context acceptor read ids.
    def get_all_variant_context_acceptor_read_ids(self):
        acceptorreadids = []
        for varcon in self.variant_contexts.values():
            acceptorreadids.extend(varcon.get_acceptor_read_ids())
        return acceptorreadids

    # Returns the variant context donor read ids.
    def get_all_variant_context_donor_read_ids(self):
        donorreadids = []
        for varcon in self.variant_contexts.values():
            donorreadids.extend(varcon.get_donor_read_ids())
        return donorreadids

    # Returns the variant context field data
    def get_variant_context_fields(self):
        return self.varcon_fields

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA=============================
    def get_variant_context_areads(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_reads()
        return []

    # Returns the donor reads of a specified variant context
    def get_variant_context_dreads(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_reads()
        return []

    # Returns the acceptor context reads
    def get_acceptor_context_reads(self, contextid):
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_acceptor_context():
                return self.variant_contexts[contextid].get_acceptor_context_reads()
        return []

    # Returns the reads of the donor context
    def get_donor_context_reads(self, contextid):
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_donor_context():
                return self.variant_contexts[contextid].get_donor_context_reads()
        return []

    # ===BASIC VARIANTCONTEXTFILE METHODS======================================
    # Reads a provided variant context file and saves data according to set filters.
    def read_variant_context_file(self, fileloc, samplefilter=None,
                                  IDfilter=None, chromfilter=None):
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

            IDpass = self.passes_filter(record[0], IDfilter)
            samplepass = self.passes_filter(record[1], samplefilter)
            chrompass = self.passes_filter(record[2], chromfilter)
            if not (IDpass and samplepass and chrompass):
                continue

            acceptor_reads = [ReadIdObject(readid) for readid in record[11].split(";")]
            donor_reads = [ReadIdObject(readid) for readid in record[12].split(";")]
            new_varcon = VariantContext(*record[:6], acceptor_reads, donor_reads)
            self.variant_contexts[record[0]] = new_varcon

    # Reads an acceptor context file
    def read_acceptor_context_file(self, accconfileloc, samplefilter=None, contextfilter=None, chromfilter=None):
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

    # Returns whether something is in the filter or not.
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
    # Returns the list of variant context identifiers from both variant context files.
    def get_variant_contexts_union(self, other_varcon_file):
        """Returns all variant context identifiers in both VariantContextFiles.

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

    # Returns the list of variant context identifiers present in both variant context files.
    def get_variant_contexts_intersect(self, other_varcon_file):
        """

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

    # Returns the list of variant context identifiers in this file but not present in the other variant context file.
    def get_variant_contexts_difference(self, other_varcon_file):
        """Determines and return the variant context identifiers not present in the other variant context file.

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

    # Returns the list of variant context identifiers only present in one of the variant context file but not the other.
    def get_variant_contexts_symmetric_difference(self, other_varcon_file):
        """Determines the variant context identifiers present in either one or the other.

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

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA BASED ON FILTERS===============
    # Returns a list/hashmap of VariantContextObjects.
    def get_variant_contexts2(self, aslist=False, varconfilter=None,
                              samplefilter=None, chromfilter=None):
        """

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
        list or dict
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
    # Main method that returns whether a variant (SNP or indel).
    def variant_is_in_context(self, varianttype, searchchrom,
                              searchstart, searchstop):
        """

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

        """
        if varianttype == "snp":
            return self.snp_variant_is_in_context(searchchrom, searchstart)
        if varianttype == "indel":
            return self.indel_variant_is_in_context(searchchrom, searchstart,
                                                    searchstop)
        return None

    # Determines whether an SNP variant is located in an already existing variant context.
    def snp_variant_is_in_context(self, varchrom, vcfvarpos):
        for varcon in self.variant_contexts.values():
            if varchrom == varcon.get_variant_context_chrom():
                if (vcfvarpos >= varcon.get_variant_context_start()
                   and vcfvarpos <= varcon.get_variant_context_end()):
                    return True
        return False

    # Determines whether an indel variant is located within an existing variant context (indelLeftPos and indelRightPos
    # can be the used search window).
    def indel_variant_is_in_context(self, indelchrom, indelleftpos, indelrightpos):
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

    # Checks whether a context is located in an already existing context
    def context_collision(self, context_arr):
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
    # Sets a variant context object to a specified context id.
    def set_variant_context(self, varconid, varcontext):
        self.variant_contexts[varconid] = varcontext

    # Adds a variant context object
    def add_variant_context(self, varconid, sampleid,
                            varconchrom, varconorigin,
                            varconstart, varconend,
                            varcon_areads, varcon_dreads,
                            acceptor_context=None, donor_context=None):
        varcon_obj = VariantContext(varconid, sampleid,
                                    varconchrom, varconorigin,
                                    varconstart, varconend,
                                    varcon_areads, varcon_dreads,
                                    acceptor_context, donor_context)
        self.variant_contexts[varconid] = varcon_obj

    # Adds an existing variant context to the variant context file
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

    # Adds an acceptor context object to a variant context.
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

    # Add a newly created acceptor context to an existing variant context.
    def add_acceptor_context(self, contextid, sampleid,
                             contextchrom, contextorigin,
                             contextstart, contextend,
                             acceptorreads):
        """Creates and adds an acceptor context to a variant context.

        Parameters
        ----------
        contextid : str
            Context identifier
        sampleid: str
            Sample name/identifier
        contextchrom: str
            Chromosome name the context is located on
        contextorigin: int
            Variant genomic position that created the context
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

    # Sets the donor context object of a variant context
    def set_donor_context(self, varconid, donor_context):
        """Adds an already existing context

        Parameters
        ----------
        varconid : str
            Context identifier
        donor_context : OverlapContext
            Donor context to add
        """
        if varconid in self.variant_contexts:
            self.variant_contexts[varconid].set_donor_context(donor_context)

    # Add a newly created donor context to an existing variant context.
    def add_donor_context(self, contextid, sampleid,
                          contextchrom, contextorigin,
                          contextstart, contextend,
                          donorreads):
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].add_donor_context(
                    contextid, sampleid,
                    contextchrom, contextorigin,
                    contextstart, contextend,
                    donorreads)

    # ===METHODS TO ADD UNMAPPED MATE IDS TO ACCEPTOR, DONOR AND VARIANT CONTEXT=======
    # Sets the unmapped mate read id for a specified acceptor context.
    def set_acceptor_context_unmapped_mate_ids(self, contextid, mateids):
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_acceptor_context_unmapped_mates(mateids)

    # Sets the unmapped mate read id for a specified donor context.
    def set_donor_context_unmapped_mate_ids(self, contextid, mateids):
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_donor_context_unmapped_mates(mateids)

    # Sets the acceptor unmapped mate ids for a specified variant context.
    def set_unmapped_acceptor_mate_ids(self, contextid, mateids):
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_acceptor_mate_ids(mateids)

    # Sets the donor unmapped mate ids for a specified variant context.
    def set_unmapped_donor_mate_ids(self, contextid, mateids):
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_donor_mate_ids(mateids)

    # ===METHODS TO GET UNMAPPED MATE IDS OF AN ACCEPTOR, DONOR AND VARIANT CONTEXT============
    # Return the unmapped read mate identifiers of a specified acceptor context.
    def get_acceptor_context_unmapped_mate_ids(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_context_unmapped_mate_ids()
        return []

    # Returns the unmapped read mate identifiers of a specified donor context.
    def get_donor_context_unmapped_mate_ids(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_context_unmapped_mate_ids()
        return []

    # Returns the unmapped acceptor read mate identifiers of a specified variant context.
    def get_unmapped_acceptor_mate_ids(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_acceptor_read_ids()
        return []

    # Returns the unmapped donor read mate identifiers of a specified variant context.
    def get_unmapped_donor_mate_ids(self, contextid):
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_donor_read_ids()
        return []

    # ===METHODS TO OBTAIN SOME STATISTICS ABOUT ALL THE CONTEXTS==============
    # Returns the average variant context length within this variant context file.
    def get_average_variant_context_length(self):
        return statistics.mean([varcon.get_variant_context_length()
                                for varcon in self.variant_contexts.values()])

    # Returns the median variant context length within this variant context file.
    def get_median_variant_context_length(self):
        return statistics.median([varcon.get_variant_context_length()
                                  for varcon in self.variant_contexts.values()])

    # Returns the average number of variant context acceptor reads
    def get_average_variant_context_acceptor_reads(self):
        return statistics.mean([varcon.get_number_of_acceptor_reads() for varcon in self.variant_contexts.values()])

    # Returns the average number of variant context donor reads
    def get_average_variant_context_donor_reads(self):
        return statistics.mean([varcon.get_number_of_donor_reads() for varcon in self.variant_contexts.values()])

    # Returns the median number of variant context acceptor reads
    def get_median_variant_context_acceptor_reads(self):
        return statistics.median([varcon.get_number_of_acceptor_reads() for varcon in self.variant_contexts.values()])

    # Returns the median number of variant context donor reads
    def get_median_variant_context_donor_reads(self):
        return statistics.median([varcon.get_number_of_donor_reads() for varcon in self.variant_contexts.values()])

    # ===METHODS TO WRITE VARIANT CONTEXT DATA TO A FILE=======================
    # Writes the variant context data to an output file.
    def write_variant_context_file(self, outfileloc, samplefilter=None,
                                   varconfilter=None, chromfilter=None):
        try:
            with open(outfileloc, "w") as varcon_outfile:
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

    # Writes the donor contexts used to construct the variant contexts.
    def write_acceptor_context_file(self, outfileloc, samplefilter=None,
                                    contextfilter=None, chromfilter=None):
        try:
            with open(outfileloc, "w") as varcon_oufile:
                varcon_oufile.write("#ContextId\tDonorSample\tChrom\tOrigin\t"
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
                        varcon_oufile.write(
                            varcon.get_acceptor_context().to_string() + "\n"
                                )
        except IOError as ioe:
            self.vaselogger.warning("Could not write acceptor contexts to "
                                    f"{ioe.filename}")

    # Writes the acceptor cotnexts used to construct the variant contexts.
    def write_donor_context_file(self, outfileloc, samplefilter=None,
                                 contextfilter=None, chromfilter=None):
        try:
            with open(outfileloc, "w") as varcon_outfile:
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

    # Writes some statistics about the acceptor and donor reads identified for each variant context.
    def write_variant_context_stats(self, statsoutloc):
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
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

    # Writes some statistics about the acceptor and donor reads identified for each variant context.
    def write_acceptor_context_stats(self, statsoutloc):
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
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
    def write_donor_context_stats(self, statsoutloc):
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
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

    # Writes the left and right positions to the output file.  Left pos for R1 and right pos for R2.
    def write_left_right_positions(self, typetowrite, outfileloc):
        try:
            with open(outfileloc, "w") as lrpof:
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

    # Writes the acceptor context left and right positions to the output file. Left pos for R1 and right pos for R2.
    def write_acceptor_left_right_positions(self, outfileloc):
        try:
            with open(outfileloc, "w") as lrpof:
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

    # Writes the left and right positions to the output file.  Left pos for R1 and right pos for R2.
    def write_donor_left_right_positions(self, outfileloc):
        try:
            with open(outfileloc, "w") as lrpof:
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
        except IOError as ioe:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    # Writes the identifiers of reads that have unmapped mates per
    # sample to a file.  Samples are all donors and the ?template?.
    def write_reads_with_unmapped_mate(self, typetowrite, umfileloc):
        """Writes variant context read identifiers with unmapped mates to a specified output file.

        Parameters
        ----------
        typetowrite : str
            Write variant context acceptor or donor read identifiers
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umFile:
                umFile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    if typetowrite == "acceptor":
                        umFile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_acceptor_mate_ids())
                            )
                    if typetowrite == "donor":
                        umFile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_donor_mate_ids())
                            )
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    # Writes the unmapped mate id of the acceptor context.
    def write_acceptor_unmapped_mates(self, umfileloc):
        """Writes the acceptor context read identifiers that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    acccon = varcon.get_acceptor_context()
                    umfile.write(str(acccon.get_context_id()) + "\t"
                                 + str(acccon.get_sample_id()) + "\t"
                                 + ";".join(acccon.get_unmapped_read_mate_ids()))
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    # Writes the unmapped mate id of the donor context.
    def write_donor_unmapped_mates(self, umfileloc):
        """Writes the donor context read identifiers that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umFile:
                umFile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    doncon = varcon.get_donor_context()
                    umFile.write(str(doncon.get_context_id()) + "\t"
                                 + str(doncon.get_sample_id()) + "\t"
                                 + ";".join(doncon.get_unmapped_read_mate_ids()))
        except IOError:
            self.vaselogger.warning("Could not write read identifiers of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    # Compares the current VariantContextFile to another variant context file
    def compare(self, othervarconfile, contextfilter=None):
        varcondiffs = {}
        for contextid in self.variant_contexts:
            if self.passes_filter(contextid, contextfilter):
                diffs = self.variant_contexts[contextid].compare(othervarconfile.get_variant_context(contextid))
                varcondiffs[contextid] = diffs
        return varcondiffs
