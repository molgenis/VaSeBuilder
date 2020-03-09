"""VariantContextFile object class.

The VariantContextFile class stores VariantContext objects and is used to
write variants contexts to an external file. This class also defines
methods accessing data from multiple VariantContexts.
"""

import logging
import statistics as stats
import sys

from variant_context import VariantContext
from overlap_context import OverlapContext
from read_id_object import ReadIdObject
from inclusion_filter import InclusionVariant


class VariantContextFile:
    """The VariantContextFile saves variant contexts.

    Attributes
    ----------
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
        self.template_alignment_file = ""
        self.contributing_alignment_files = []
        self.contributing_variant_files = []

        # Check whether to read a provided variant context file with set optional parameters
        if fileloc is not None:
            self.read_variant_context_file(fileloc, samplefilter, varconfilter, chromfilter)

    # ===METHODS TO GET DATA FROM THE VARIANT CONTEXT FILE=====================
    def get_variant_contexts(self, asdict=False):
        """Return all variant contexts.

        By default all variant contexts are gathered in a list and returned.
        If the asdict parameter is set to True, the variant contexts are
        returned as a dictionary.

        Returns
        -------
        list or dict
            Variant context as list by default, as dict when parameter is set to True
        """
        if asdict:
            return self.variant_contexts
        return list(self.variant_contexts.values())

    def get_variant_contexts_by_sampleid(self):
        """Return all variant contexts as dictionary grouped by sample name.

        Returns
        -------
        dict

        """
        varcons = self.get_variant_contexts()
        return {x.get_variant_context_sample(): [
            y for y in varcons
            if y.get_variant_context_sample() == x.get_variant_context_sample()
            ] for x in varcons}

    def get_number_of_contexts(self):
        """Count variant contexts.

        Returns
        -------
        int
            Number of variant contexts
        """
        return len(self.variant_contexts)

    def get_variant_context_ids(self):
        """Return all variant context identifiers.

        Returns
        -------
        list of str
            Context identifiers of all variant contexts
        """
        return list(self.variant_contexts.keys())

    def get_variant_context(self, contextid):
        """Return a variant context using its ID, if the context exists.

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
        """Check if a variant context ID exists.

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
        """Return the acceptor context of the specified variant context.

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
        """Return the donor context of the specified variant context.

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
        """Return all acceptor reads from all variant contexts.

        Returns
        -------
        list of pysam.AlignedSegment
            Acceptor reads of all variant contexts
        """
        acceptorreads = []
        for varcon in self.variant_contexts.values():
            acceptorreads.extend(varcon.get_acceptor_reads())
        return acceptorreads

    def get_all_variant_context_donor_reads(self):
        """Return all donor reads from all variant contexts as string tuples.

        Returns
        -------
        list of tuples
            Variant context donor reads
        """
        dbrs = []
        donorreads = []
        for varcon in self.variant_contexts.values():
            dbrs.extend(varcon.get_donor_reads())
        for dbr in dbrs:
            readpn = "2"
            if dbr.is_read1:
                readpn = "1"
            donorreads.append((dbr.query_name,
                               readpn,
                               dbr.query_sequence,
                               "".join([chr(x+33) for x in dbr.query_qualities])))
        return list(set(donorreads))

    def get_all_variant_context_donor_reads_2(self):
        """Return all donor reads from all variant contexts.

        Returns
        -------
        donor_reads : list of pysam.AlignedSegment objects
            List of all reads in all variant contexts
        """
        donor_reads = []
        for varcon in self.variant_contexts.values():
            donor_reads.extend(varcon.get_donor_reads())
        return donor_reads

    def get_all_variant_context_acceptor_read_ids(self):
        """Return all variant context acceptor read IDs.

        Returns
        -------
        acceptorreadids : list of str
            Variant context acceptor read IDs
        """
        acceptorreadids = []
        for varcon in self.variant_contexts.values():
            acceptorreadids.extend(varcon.get_acceptor_read_ids())
        return acceptorreadids

    def get_all_variant_context_donor_read_ids(self):
        """Return all variant context donor read IDs.

        Returns
        -------
        donorreadids : list of str
            Variant context donor read IDs
        """
        donorreadids = []
        for varcon in self.variant_contexts.values():
            donorreadids.extend(varcon.get_donor_read_ids())
        return donorreadids

    def get_all_variant_context_variant_records(self):
        """Return all variant records from all variant contexts, sorted."""
        variants = []
        for varcon in self.variant_contexts.values():
            variants.extend(varcon.variants)
        variants.sort(key=lambda x: x.pos)
        chr_sort_order = [str(x) for x in range(1, 23)] + ["X", "Y", "MT"]
        variants.sort(key=lambda x: chr_sort_order.index(x.chrom))
        return variants

    def get_variant_context_fields(self):
        """Return the number representations of the variant context data.

        Returns
        -------
        self.varcon_fields : dict
            Number representation of each variant context field
        """
        return self.varcon_fields

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA================================
    def get_variant_context_areads(self, contextid):
        """Return the variant context acceptor reads of a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of pysam.AlignedSegment
            Variant context acceptor reads, empty list of context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_reads()
        return []

    def get_variant_context_dreads(self, contextid):
        """Return the variant context donor reads of a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of pysam.AlignedSegment
            Variant context donor reads, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_reads()
        return []

    def get_acceptor_context_reads(self, contextid):
        """Return the reads of a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Acceptor context identifier

        Returns
        -------
        list of pysam.AlignedSegment
            Acceptor context reads
        """
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_acceptor_context():
                return self.variant_contexts[contextid].get_acceptor_context_reads()
        return []

    def get_donor_context_reads(self, contextid):
        """Return the reads of a specified donor context.

        Parameters
        ----------
        contextid : str
            Donor context identifier

        Returns
        -------
        list of pysam.AlignedSegment
            Donor context reads, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            if self.variant_contexts[contextid].has_donor_context():
                return self.variant_contexts[contextid].get_donor_context_reads()
        return []

    # ===BASIC VARIANTCONTEXTFILE METHODS======================================
    def read_variant_context_file(self, fileloc, samplefilter=None,
                                  idfilter=None, chromfilter=None):
        """Read a provided variant context file and save the data.

        Filter can be set for reading the variant context file. The sample
        filter can be used to specify which samples to save. Samples not in
        the samplefilter will be skipped. Similarly filter for variant
        contexts and chromosome names can be set.

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
            sys.exit()
        varcon_records = [x.strip().split("\t")
                          for x in varcon_records if not x.startswith("#")]

        for record in varcon_records:
            if record[0] in self.variant_contexts:
                continue

            id_pass = self.passes_filter(record[0], idfilter)
            samplepass = self.passes_filter(record[1], samplefilter)
            chrompass = self.passes_filter(record[2], chromfilter)
            if not (id_pass and samplepass and chrompass):
                continue

            acceptor_reads = [ReadIdObject(readid) for readid in record[11].split(";")]
            donor_reads = [ReadIdObject(readid) for readid in record[12].split(";")]
            donor_variants = [InclusionVariant(record[1], *var.split("_"))
                              for var in record[13].split(";")]
            new_varcon = VariantContext(*record[:6], acceptor_reads, donor_reads,
                                        variants=donor_variants)
            self.variant_contexts[record[0]] = new_varcon

    def read_acceptor_context_file(self, accconfileloc, samplefilter=None,
                                   contextfilter=None, chromfilter=None):
        """Add a provided acceptor context from file to already existing variant contexts.

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
                filelines = accconfile.readlines()
        except IOError:
            self.vaselogger.warning(f"Could not read acceptor context file: {accconfileloc}")

        filelines = [line.strip().split("\t") for line in filelines
                     if not line.startswith("#")]

        for fileline in filelines:
            samplepass = self.passes_filter(fileline[1], samplefilter)
            contextpass = self.passes_filter(fileline[0], contextfilter)
            chrompass = self.passes_filter(fileline[2], chromfilter)
            if not (samplepass and contextpass and chrompass):
                continue
            if not fileline[0] in self.variant_contexts:
                continue

            self.variant_contexts[fileline[0]].add_acceptor_context(
                fileline[0], fileline[1], fileline[2], int(fileline[3]),
                int(fileline[4]), int(fileline[5]), fileline[7].split(";")
                )

    # Reads a donor context file.
    def read_donor_context_file(self, donconfileloc, samplefilter=None,
                                contextfilter=None, chromfilter=None):
        """Add a provided donor context from file to already existing variant contexts.

        Parameters
        ----------
        donconfileloc : str
            Path to the donor contexts file
        samplefilter : list of str
            Sample names/identifiers to include
        contextfilter : list of str
            Donor contexts to include
        chromfilter : list of str
            Chromosome names to include
        """
        try:
            with open(donconfileloc, "r") as donconfile:
                filelines = donconfile.readlines()
        except IOError:
            self.vaselogger.warning(f"Could not read acceptor context file: {donconfileloc}")

        filelines = [line.strip().split("\t") for line in filelines
                     if not line.startswith("#")]

        for fileline in filelines:
            samplepass = self.passes_filter(fileline[1], samplefilter)
            contextpass = self.passes_filter(fileline[0], contextfilter)
            chrompass = self.passes_filter(fileline[2], chromfilter)
            if not (samplepass and contextpass and chrompass):
                continue
            if not fileline[0] in self.variant_contexts:
                continue

            self.variant_contexts[fileline[0]].add_donor_context(
                fileline[0], fileline[1], fileline[2], int(fileline[3]),
                int(fileline[4]), int(fileline[5]), fileline[7].split(";")
                )

    @staticmethod
    def passes_filter(valtocheck, filterlist):
        """Test if a provided value passes a provided inclusion filter.

        Filters are expected to be inclusion filters, denoting a list of values
        that should be used/included.

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
        """Return all variant context identifiers from both variant context files.

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
        """Return variant context identifiers present in both variant context files.

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
        """Return the variant context identifiers not present in the other variant context file.

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
        """Return 'xor' of variant contexts between two variant context files.

        THe normal difference only returns the context identifiers in A but not
        in B. The symmetric difference returns the context identifiers in A but
        not in B as well as context identifiers in B but not in A.

        Parameters
        ----------
        other_varcon_file : VariantContextFile

        Returns
        -------
        list of str
            Variant context read IDs
        """
        own_varcon_ids = self.get_variant_context_ids()
        other_varcon_ids = other_varcon_file.get_variant_context_ids()
        return list(set(own_varcon_ids) ^ set(other_varcon_ids))

    # ===METHODS TO OBTAIN VARIANT CONTEXT DATA BASED ON FILTERS===============.
    def get_variant_contexts2(self, aslist=False, varconfilter=None,
                              samplefilter=None, chromfilter=None):
        """Return all variant contexts in the VariantContextFile.

        Inclusion filters can be set to subset the returned variant contexts.
        Filters include sample names/identifiers, variant context identifiers,
        and chromosome names. Variant contexts can be returned as list or
        dictionary.

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
            return [x for x in self.variant_contexts.values()
                    if (self.passes_filter(x.get_variant_context_id(),
                                           varconfilter)
                        and self.passes_filter(x.get_variant_context_sample(),
                                               samplefilter)
                        and self.passes_filter(x.get_variant_context_chrom(),
                                               chromfilter))]
        return {k: v for k, v in self.variant_contexts
                if (self.passes_filter(k, varconfilter)
                    and self.passes_filter(v.get_variant_context_sample(),
                                           samplefilter)
                    and self.passes_filter(v.get_variant_context_chrom(),
                                           chromfilter))}

    # ====METHODS TO ASSESS WHETHER A VARIANT IS IN AN EXISTING CONTEXT========
    def variant_is_in_context(self, varianttype, searchchrom,
                              searchstart, searchstop):
        """Check if a variant is located in an already existing context.

        If the variant type is set to 'snp' the method will check and return
        whether the SNP overlaps with a variant context. If the variant is set
        to 'indel' the method will check and return whether the area from start
        to end of the indel overlaps with a variant context.

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
        VariantContext or None
            The variant overlapping with the variant, None if no overlap
        """
        if varianttype == "snp":
            return self.snp_variant_is_in_context(searchchrom, searchstart)
        if varianttype == "indel":
            return self.indel_variant_is_in_context(searchchrom, searchstart,
                                                    searchstop)
        return None

    def snp_variant_is_in_context(self, varchrom, vcfvarpos):
        """Check if a SNP is located in an already existing variant context.

        For each variant context it is first checked whether the chromosome
        name of the context and SNP are the same. If so, it is then checked
        whether the SNP genomic position is equal or larger than the context
        start and smaller or equal than the context end.

        Parameters
        ----------
        varchrom : str
            Chromosome name the SNP is located on
        vcfvarpos : int
            Genomic position of the SNP
        Returns
        -------
        VariantContext or None
            Variant context overlapping with the variant, None if no context overlaps
        """
        for varcon in self.variant_contexts.values():
            if varchrom == varcon.get_variant_context_chrom():
                start = varcon.get_variant_context_start()
                stop = varcon.get_variant_context_end()
                if start <= vcfvarpos <= stop:
                    return varcon
        return None

    def indel_variant_is_in_context(self, indelchrom, indelleftpos, indelrightpos):
        """Check if an indel is located in an already existing variant context.

        For each variant context, first checks whether the chromosome name of
        the context and indel are the same. If so, it is then checked whether
        the indel range from indel start to end overlaps with the context.

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
        VariantContext or None
            Variant context overlapping with the variant, None if no overlap
        """
        for varcon in self.variant_contexts.values():
            if indelchrom == varcon.get_variant_context_chrom():
                if (indelleftpos <= varcon.get_variant_context_start()
                        and indelrightpos >= varcon.get_variant_context_start()):
                    return varcon
                if (indelleftpos <= varcon.get_variant_context_end()
                        and indelrightpos >= varcon.get_variant_context_end()):
                    return varcon
                if (indelleftpos >= varcon.get_variant_context_start()
                        and indelrightpos <= varcon.get_variant_context_end()):
                    return varcon
        return None

    def context_collision(self, context_arr):
        """Check if a potential context overlaps with an existing context.

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
        for varcon in self.variant_contexts.values():
            if varcon.get_variant_context_chrom() == context_arr[0]:
                if (varcon.get_variant_context_start() <= context_arr[3]
                        and context_arr[2] <= varcon.get_variant_context_end()):
                    return True
        return False

    def context_collision_v2(self, context_arr):
        """Check if a potential context overlaps with an existing context.

        Parameters
        ----------
        context_arr: list
            Essential context data (chrom, variant pos, start, end)

        Returns
        -------
        VariantContext or None
            Variant context in which the overlap occurs, None if there is no overlap
        """
        if f"{context_arr[0]}_{context_arr[1]}" in self.variant_contexts:
            return self.variant_contexts[f"{context_arr[0]}_{context_arr[1]}"]
        for varcon in self.variant_contexts.values():
            if varcon.get_variant_context_chrom() == context_arr[0]:
                if (varcon.get_variant_context_start() <= context_arr[3]
                        and context_arr[2] <= varcon.get_variant_context_end()):
                    return varcon
        return None

    # ===METHODS TO ADD DATA/VARIANT CONTEXTS TO THE VARIANT CONTEXT FILE======
    def set_variant_context(self, varconid, varcontext):
        """Set a provided variant context with the provided context identifier.

        Will overwrite the previous value for the provided context identifier.

        Parameters
        ----------
        varconid : str
            Identifier of the variant context to add
        varcontext: VariantContext
            The VariantContext to set
        """
        self.variant_contexts[varconid] = varcontext

    def set_variant_context_donor_reads(self, varconid, donor_reads):
        """Set the variant context donor reads for a specified variant context.

        Parameters
        ----------
        varconid : str
            Variant context identifier to set donor reads for
        donor_reads : list of DonorBamRead
            List of donor reads to set
        """
        if varconid in self.variant_contexts:
            self.variant_contexts[varconid].set_donor_reads(donor_reads)

    def add_variant_context(self, varconid, sampleid,
                            varconchrom, varconorigin,
                            varconstart, varconend,
                            varcon_areads, varcon_dreads,
                            acceptor_context=None, donor_context=None):
        """Construct and add a VariantContext from the provided parameters.

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
            Variant context donor reads
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

        The variant context to add should be a valid VariantContext object.
        It will be added to the variant context file under the context
        identifier.

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
        """Add an existing acceptor context to a variant context.

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
        """Construct and add an acceptor context to a specified variant context.

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
                acceptorreads
                )

    def set_donor_context(self, varconid, donor_context):
        """Add an already existing donor context to a specified variant context.

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
        """Construct and add a donor context to a specified variant context.

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
                donorreads
                )

    # ===METHODS TO ADD UNMAPPED MATE IDS TO ACCEPTOR, DONOR AND VARIANT CONTEXT=======
    def set_acceptor_context_unmapped_mate_ids(self, contextid, mateids):
        """Set the read IDs that have unmapped read mates for a specified acceptor context.

        Parameters
        ----------
        contextid : str
            Acceptor context to set unmapped read IDs for
        mateids : lit of str
            read IDs with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_acceptor_context_unmapped_mates(mateids)

    def set_donor_context_unmapped_mate_ids(self, contextid, mateids):
        """Set the read IDs that have unmapped read mates for a specified donor context.

        Parameters
        ----------
        contextid : str
            Donor context to set unmapped read IDs for
        mateids : list of str
            Donor context read IDs with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_donor_context_unmapped_mates(mateids)

    def set_unmapped_acceptor_mate_ids(self, contextid, mateids):
        """Set the acceptor read IDs that have unmapped mates for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier to set unmapped acceptor read identifiers for
        mateids : list of str
            Variant context acceptor read IDs with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_acceptor_mate_ids(mateids)

    def set_unmapped_donor_mate_ids(self, contextid, mateids):
        """Set the donor read IDs that have unmapped mates for a specified variant context.

        Parameters
        ----------
        contextid : str
            Variant context identifier to set unmapped donor read IDs for
        mateids: list of str
            Variant context donor read IDs with unmapped mates
        """
        if contextid in self.variant_contexts:
            self.variant_contexts[contextid].set_unmapped_donor_mate_ids(mateids)

    # ===METHODS TO GET UNMAPPED MATE IDS OF AN ACCEPTOR, DONOR AND VARIANT CONTEXT============
    def get_acceptor_context_unmapped_mate_ids(self, contextid):
        """Return the acceptor context read IDs that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Acceptor context identifier

        Returns
        -------
        list of str
            Acceptor context read IDs, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_acceptor_context_unmapped_mate_ids()
        return []

    def get_donor_context_unmapped_mate_ids(self, contextid):
        """Return the donor context read IDs that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Donor context identifier

        Returns
        -------
        list of str
            Donor contexts read IDs, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_donor_context_unmapped_mate_ids()
        return []

    def get_unmapped_acceptor_mate_ids(self, contextid):
        """Return the variant context acceptor read IDs that have unmapped read mates.

        Parameters
        ----------
        contextid : str
            Variant context identifier

        Returns
        -------
        list of str
            Variant context acceptor read identifiers, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_acceptor_read_ids()
        return []

    def get_unmapped_donor_mate_ids(self, contextid):
        """Return the variant context donor read IDs that have unmapped read mates.

        Parameters
        ----------
        contextid : str

        Returns
        -------
        list of str
            Variant context donor read IDs, empty list if context does not exist
        """
        if contextid in self.variant_contexts:
            return self.variant_contexts[contextid].get_unmapped_donor_read_ids()
        return []

    # ===METHODS TO OBTAIN SOME STATISTICS ABOUT ALL THE CONTEXTS==============
    def get_average_variant_context_length(self):
        """Calculate the mean variant context length.

        The mean length is calculated over all variant contexts in the variant
        context file. Acceptor and donor contexts associated with the variant
        contexts are not used.

        Returns
        -------
        float
            Average variant context length
        """
        return stats.mean([varcon.get_variant_context_length()
                           for varcon in self.variant_contexts.values()])

    def get_median_variant_context_length(self):
        """Calculate the median variant context length.

        The median length is calculated over all variant contexts in the
        variant context file. Acceptor and donor contexts associated with the
        variant contexts are not used.

        Returns
        -------
        float
            Median variant context length
        """
        return stats.median([varcon.get_variant_context_length()
                             for varcon in self.variant_contexts.values()])

    def get_average_variant_context_acceptor_reads(self):
        """Calculate the average number of variant context acceptor reads.

        The mean number of acceptor reads is calculated over all variant contexts.

        Returns
        -------
        float
            Mean number of variant context acceptor reads
        """
        return stats.mean([varcon.get_number_of_acceptor_reads()
                           for varcon in self.variant_contexts.values()])

    def get_average_variant_context_donor_reads(self):
        """Calculate the average number of variant context donor reads.

        The mean number of donor reads is calculated over all variant contexts.

        Returns
        -------
        float
            Mean number of variant context donor reads
        """
        return [varcon.get_number_of_donor_reads()
                for varcon in self.variant_contexts.values()]

    def get_median_variant_context_acceptor_reads(self):
        """Calculate the median variant context acceptor reads.

        Returns
        -------
        float
            Median number of variant context acceptor reads
        """
        return stats.median([varcon.get_number_of_acceptor_reads()
                             for varcon in self.variant_contexts.values()])

    def get_median_variant_context_donor_reads(self):
        """Calculate the median number of variant context donor reads.

        Returns
        -------
        float
            Median number of variant context donor reads
        """
        return stats.median([varcon.get_number_of_donor_reads()
                             for varcon in self.variant_contexts.values()])

    # ===METHODS TO WRITE VARIANT CONTEXT DATA TO A FILE=======================
    def write_variant_context_file(self, outfileloc, vbuuid, samplefilter=None,
                                   varconfilter=None, chromfilter=None):
        """Write the data saved in the variant contexts to an output file.

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
                                     "DonorReadIds\tDonorVariants\n")
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
        """Write the acceptor contexts associated with variants contexts to an output file.

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
                varcon_outfile.write(f"#VBUUID: {vbuuid}\n")
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
        """Write the donor contexts associated with variant context to an output file.

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
                varcon_outfile.write(f"#VBUUID: {vbuuid}\n")
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
        """Write basic statistics of the variant contexts to a specified output file.

        Includes means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc : str
            Path and name to write variant context statistics output file to
        vbuuid : str
            VaSeBuilder identifier
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"#VBUUID: {vbuuid}\n")
                varcon_statsfile.write("#ContextId\tAvg_ALen\tAvg_DLen\t"
                                       "Med_ALen\tMed_DLen\tAvg_AQual\t"
                                       "Avg_DQual\tMed_AQual\tMed_DQual\t"
                                       "Avg_AMapQ\tAvg_DMapQ\tMed_AMapQ\t"
                                       "Med_DMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(varcon.to_statistics_string() + "\n")
        except IOError:
            self.vaselogger.critical("Could not write variant context "
                                     f"statistics to {statsoutloc}")

    def write_acceptor_context_stats(self, statsoutloc, vbuuid):
        """Write basic statistics of the acceptor contexts to a specified output file.

        Includes means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc: str
            Path and name to write the acceptor context statistics file to
        vbuuid : str
            VaSeBuilder identifier
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"#VBUUID: {vbuuid}\n")
                varcon_statsfile.write("#ContextId\tAvg_ReadLen\tMed_ReadLen\t"
                                       "Avg_ReadQual\tMed_ReadQual\t"
                                       "Avg_ReadMapQ\tMed_ReadMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(
                        varcon.get_acceptor_context().to_statistics_string()
                        + "\n"
                        )
        except IOError:
            self.vaselogger.critical("Could not write acceptor context "
                                     f"statistics to {statsoutloc}")

    # Writes some statistics about the acceptor and donor reads identified for each variant context.
    def write_donor_context_stats(self, statsoutloc, vbuuid):
        """Write basic statistics of the donor contexts to a specified output file.

        Includes means and medians of read lengths, Q-Score and mapping quality.

        Parameters
        ----------
        statsoutloc : str
            Path and name to write the donor context statistics to
        """
        try:
            with open(statsoutloc, "w") as varcon_statsfile:
                varcon_statsfile.write(f"#VBUUID: {vbuuid}\n")
                varcon_statsfile.write("#ContextId\tAvg_ReadLen\tMed_ReadLen\t"
                                       "Avg_ReadQual\tMed_ReadQual\t"
                                       "Avg_ReadMapQ\tMed_ReadMapQ\n")
                for varcon in self.variant_contexts.values():
                    varcon_statsfile.write(
                        varcon.get_donor_context().to_statistics_string()
                        + "\n"
                        )
        except IOError:
            self.vaselogger.critical("Coud not write donor context statistics "
                                     f"to {statsoutloc}")

    def write_left_right_positions(self, typetowrite, outfileloc, vbuuid):
        """Write variant context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic positions of R1 reads and rightmost genomic
        positions for R2 reads are written to file.

        Parameters
        ----------
        typetowrite : str
            Type (acceptor/donor) of file to write
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"#VBUUID: {vbuuid}")
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
        except IOError:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    def write_acceptor_left_right_positions(self, outfileloc, vbuuid):
        """Write acceptor context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic position of R1 reads and rightmost genomic
        positions for R2 reads are written to file.

        Parameters
        ----------
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"#VBUUID: {vbuuid}\n")
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
        except IOError:
            self.vaselogger.warning("Could not write read left positions to "
                                    f"output file {outfileloc}")

    def write_donor_left_right_positions(self, outfileloc, vbuuid):
        """Write donor context leftmost and rightmost genomic read positions to an output file.

        The leftmost genomic position of R1 reads and rightmost genomic
        positions for R2 reads are written to file.

        Parameters
        ----------
        outfileloc : str
            Path and name to write the output file to
        """
        try:
            with open(outfileloc, "w") as lrpof:
                lrpof.write(f"#VBUUID: {vbuuid}\n")
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
        """Write variant context read IDs with unmapped mates to a specified output file.

        Parameters
        ----------
        typetowrite : str
            Write variant context acceptor or donor read IDs
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"#VBUUID: {vbuuid}\n")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    if typetowrite == "acceptor":
                        umfile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_acceptor_mate_ids())
                            + "\n")
                    if typetowrite == "donor":
                        umfile.write(
                            varcon.get_variant_context_id() + "\t"
                            + str(varcon.get_variant_context_sample()) + "\t"
                            + ";".join(varcon.get_unmapped_donor_mate_ids())
                            + "\n")
        except IOError:
            self.vaselogger.warning("Could not write read IDs of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def write_acceptor_unmapped_mates(self, umfileloc, vbuuid):
        """Write the acceptor context read IDs that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"#VBUUID: {vbuuid}\n")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    acccon = varcon.get_acceptor_context()
                    umfile.write(str(acccon.get_context_id()) + "\t"
                                 + str(acccon.get_sample_id()) + "\t"
                                 + ";".join(acccon.get_unmapped_read_mate_ids()) + "\n")
        except IOError:
            self.vaselogger.warning("Could not write read IDs of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def write_donor_unmapped_mates(self, umfileloc, vbuuid):
        """Write the donor context read IDs that have unmapped mates to an output file.

        Parameters
        ----------
        umfileloc : str
            Path to write the output file to
        """
        try:
            with open(umfileloc, "w") as umfile:
                umfile.write(f"#VBUUID: {vbuuid}\n")
                umfile.write("#ContextId\tSampleId\tReadIds\n")
                for varcon in self.variant_contexts.values():
                    doncon = varcon.get_donor_context()
                    umfile.write(str(doncon.get_context_id()) + "\t"
                                 + str(doncon.get_sample_id()) + "\t"
                                 + ";".join(doncon.get_unmapped_read_mate_ids())+"\n")
        except IOError:
            self.vaselogger.warning("Could not write read IDs of "
                                    "reads with unmapped mates to "
                                    f"{umfileloc}")

    def compare(self, othervarconfile, contextfilter=None):
        """Compare the current variant context file to another provided variant context file.

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
                diffs = self.variant_contexts[contextid].compare(
                    othervarconfile.get_variant_context(contextid)
                    )
                varcondiffs[contextid] = diffs
        return varcondiffs

    def add_variant_context_file(self, variantcontextfile):
        """Add another VariantContextFile to the existing one.

        Variant contexts from the second file overlapping with variant contexts
        from the current file will not be added.

        Parameters
        ----------
        variantcontextfile : VariantContextFile
            Variant context file to add to the current
        """
        for contextid, varcon in variantcontextfile.get_variant_contexts(True).items():
            if contextid not in self.variant_contexts:
                if not self.context_collision(varcon.get_context()):
                    self.variant_contexts[contextid] = varcon

    @staticmethod
    def merge_context_windows(context1, context2):
        """Merge two context windows.

        Parameters
        ----------
        context1 : list of str and int
            Essential data of the first context
        context2 : list of str and int
            Essential data of the second context

        Returns
        -------
        new_context_window : [chrom, origin, start, stop]
            Essential data of the combined context
        """
        new_start_pos = min([context1[2], context2[2]])
        new_end_pos = max([context1[3], context2[3]])
        new_context_window = [context1[0], context1[1], new_start_pos, new_end_pos]
        return new_context_window

    @staticmethod
    def merge_variant_context_reads(variantreads):
        """Merge the reads of two contexts.

        Parameters
        ----------
        variantreads : list of pysam.AlignedSegment
            List of reads (provided as context1_reads + context2_reads)

        Returns
        -------
        unique_variantreads : list of pysam.AlignedSegment objects
        """
        unique_variantreads = []
        checklist = []
        for varread in variantreads:
            readpn = "2"
            if varread.is_read1:
                readpn = "1"
            id_pair = (varread.query_name, readpn)
            if id_pair not in checklist:
                unique_variantreads.append(varread)
                checklist.append(id_pair)
        return unique_variantreads

    def merge_variant_contexts(self, variantcontext1, variantcontext2):
        """Merge two variant contexts and add the new context to the variant context file.

        Parameters
        ----------
        variantcontext1 : VariantContext
            First variant context to merge
        variantcontext2 : VariantContext
            Second variant context to merge
        """
        combined_acceptor_context = self.merge_overlap_contexts(
            variantcontext1.get_acceptor_context(),
            variantcontext2.get_acceptor_context()
            )
        combined_donor_context = self.merge_overlap_contexts(
            variantcontext1.get_donor_context(),
            variantcontext2.get_donor_context()
            )

        # Obtain a list of acceptor and donor reads from both variant contexts
        vareads = variantcontext1.get_acceptor_reads() + variantcontext2.get_acceptor_reads()
        vdreads = variantcontext1.get_donor_reads() + variantcontext2.get_donor_reads()

        # Combine the two variant contexts by determining the new context
        # window and acceptor and donor reads.
        combined_window = self.merge_context_windows(variantcontext1.get_context(),
                                                     variantcontext2.get_context())
        combined_vareads = self.merge_variant_context_reads(vareads)
        combined_vdreads = self.merge_variant_context_reads(vdreads)

        # Set the new combined variant context
        combined_varcon = VariantContext(
            variantcontext1.get_variant_context_id(),
            variantcontext1.get_variant_context_sample(), *combined_window,
            combined_vareads, combined_vdreads, combined_acceptor_context,
            combined_donor_context
            )
        self.set_variant_context(variantcontext1.get_variant_context_id(), combined_varcon)

    def merge_overlap_contexts(self, overlapcontext1, overlapcontext2):
        """Merge two acceptor/donor contexts and returns the combined context.

        Parameters
        ----------
        overlapcontext1 : OverlapContext
            First acceptor/donor context to merge
        overlapcontext2 : OverlapContext
            Second acceptor/donor context to merge

        Returns
        -------
        combined_accdon_context : OverlapContext
            New combined acceptor/donor context
        """
        adreads = (overlapcontext1.get_context_bam_reads()
                   + overlapcontext2.get_context_bam_reads())
        combined_window = self.merge_context_windows(overlapcontext1.get_context(),
                                                     overlapcontext2.get_context())
        combined_adreads = self.merge_variant_context_reads(adreads)
        combined_accdon_context = OverlapContext(overlapcontext1.get_context_id(),
                                                 overlapcontext1.get_sample_id(),
                                                 *combined_window, combined_adreads)
        return combined_accdon_context

    def remove_variant_context(self, contextid):
        """Remove a specified variant context.

        Parameters
        ----------
        contextid : str
            Identifier of the context to remove
        """
        if contextid in self.variant_contexts:
            del self.variant_contexts[contextid]

    def get_template_alignment_file(self):
        """Return the associated acceptor alignment file path.

        Returns
        -------
        self.template_alignment_file : str
            Path to acceptor alignment file associated with the variant context file
        """
        return self.template_alignment_file

    def set_template_alignment_file(self, alignment_file):
        """Set the acceptor alignment file associated with the variant context file.

        Parameters
        ----------
        alignment_file : str
            Path to acceptor alignment file
        """
        self.template_alignment_file = alignment_file

    def get_donor_alignment_files(self):
        """Return the associated donor alignment file paths.

        Returns
        -------
        self.contributing_alignment_files : list of str
            Donor alignment files that contributed
        """
        return self.contributing_alignment_files

    def add_donor_alignment_file(self, donor_alnfile):
        """Add an associated contributing donor alignment file.

        Parameters
        ----------
        donor_alnfile : str
            Location to contributing donor alignment file
        """
        self.contributing_alignment_files.append(donor_alnfile)

    def set_donor_alignment_files(self, donor_alnfiles):
        """Set the list of of contributing donor alignment files.

        Parameters
        ----------
        donor_alnfiles : list of str
            Locations to contributing donor alignment files
        """
        self.contributing_alignment_files = donor_alnfiles

    def get_donor_variant_files(self):
        """Return the associated donor variant files.

        Returns
        -------
        self.contributing_variant_files : list of str
            Donor variant files that contributed
        """
        return self.contributing_variant_files

    def add_donor_variant_file(self, donor_varfile):
        """Add a donor variant files.

        Parameters
        ----------
        donor_varfile : str
            Location to contributing donor variant file
        """
        self.contributing_variant_files.append(donor_varfile)

    def set_donor_variant_files(self, donor_varfiles):
        """Set the associated donor variant files.

        Parameters
        ----------
        donor_varfiles : list of str
            Donor variant files that contributed
        """
        self.contributing_variant_files = donor_varfiles
