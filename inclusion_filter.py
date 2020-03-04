"""Module to parse and subset a variant filter list."""

from collections import namedtuple

Filter = namedtuple("Filter", ["name", "values", "type"])


class InclusionFilter:
    def __init__(self, filter_file=None, subset_filters=None, priorities=None):
        self.filter_file = filter_file
        self.subset_filters = subset_filters
        self.priorities = priorities
        self.inclusions = {}

        if self.filter_file is not None:
            self.inclusions = self.read_variant_filter_file_v2(self.filter_file,
                                                               self.subset_filters,
                                                               self.priorities)

    @classmethod
    def read_variant_filter_file_v2(cls, inclusion_file, subsets=None, priorities=None):
        """

        Parameters
        ----------
        variant_filter_loc : str
            Path to variant filter file
        priority_filter : str
            Column name to select variants
        priority_values : list of str
        selection_filter : str
            Column name to use for selecting variants
        selection_values : list of str

        Returns
        -------
        variant_filter_data : dict
        """
        try:
            with open(inclusion_file) as infile:
                headers = next(infile)
                variants = infile.readlines()
        except IOError:
            return "PROBLEM"

        headers = [x.lower() for x in headers.strip().split("\t")]
        variants = [x.strip().split("\t") for x in variants]

        subsets = cls.make_filters2(subsets, "subset", headers)
        priorities = cls.make_filters2(priorities, "priority", headers)

        if len(headers) > 5:
            variants = [InclusionVariant(*x[:5], **dict(zip(headers[5:], x[5:])))
                        for x in variants]
        else:
            variants = [InclusionVariant(*x[:5]) for x in variants]

        inclusion_filter = {}
        for variant in variants:
            if subsets is not None:
                if not cls.pass_subset2(variant, subsets):
                    continue
            if priorities is not None:
                variant.priorities = cls.get_priority_levels(variant, priorities)

            if variant.sample not in inclusion_filter:
                inclusion_filter[variant.sample] = []
            inclusion_filter[variant.sample].append(variant)
        return inclusion_filter

    @staticmethod
    def get_priority_levels(variant, priorities):
        var_priorities = []
        for priority in priorities:
            value = variant.__getattribute__(priority.name)
            if priority.name == "size":
                if priority.values[0] == "greater":
                    level = variant.size
                elif priority.values[0] == "less":
                    level = float(1/variant.size)
            elif value in priority.values:
                levels = len(priority.values)
                level = levels - priority.values.index(value)
            else:
                level = 0
            var_priorities.append((priority.name, level))
        return var_priorities

    @staticmethod
    def pass_subset(variant, subset_filters):
        for subset in subset_filters:
            if variant[subset.column] in subset.values:
                continue
            return False
        return True

    @staticmethod
    def pass_subset2(variant, subset_filters):
        for subset in subset_filters:
            if subset.name == "size":
                minsize = int(subset.values[0])
                if variant.size < minsize:
                    return False
                if len(subset.values) > 1:
                    maxsize = int(subset.values[1])
                    if variant.size > maxsize:
                        return False
                return True
            elif variant.__getattribute__(subset.name) not in subset.values:
                return False
        return True

    @classmethod
    def make_filters(cls, filter_list, filter_type, header_list):
        if filter_list is None:
            return None
        new_filter_list = []
        for filt in filter_list:
            col = cls.get_filter_header_pos(filt[0], header_list)
            if col is None:
                continue
            new_filter_list.append(Filter(filt[0], filter_type, col, filt[1]))
        if not new_filter_list:
            return None
        return new_filter_list


    @classmethod
    def make_filters2(cls, filter_list, filter_type, header_list):
        if filter_list is None:
            return None
        new_filter_list = []
        allowed = header_list + ["type", "size"]
        for filt in filter_list:
            if filt[0] not in allowed:
                continue
            new_filter_list.append(Filter(filt[0], filt[1], filter_type))
        if not new_filter_list:
            return None
        return new_filter_list


    @staticmethod
    def get_filter_header_pos(filtername, headers):
        """Check and return whether

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
        if filtername is None:
            return None
        if filtername not in headers:
            return None
        return headers.index(filtername)



class InclusionVariant:
    """Saves necessary data of a genomic variant from the inclusion filter file.

    Attributes
    ----------
    vcf_variant_chrom : str
        Chromosome name the variant is located on
    vcf_variant_start : int
        The leftmost genomic position of the variant
    vcf_variant_ref : str
        The reference allele of the variant
    vcf_variant_alts : tuple
        The alternative allele(s) of the variant
    vcf_variant_filter : str
        Filter applied to the variant (i.e. PASS)
    vcf_variant_type : str
        Type of the variant (SNP/Indel/etc)
    """

    def __init__(self, sample, chrom, pos, ref, alts, **kwargs):
        self.sample = sample
        self.chrom = chrom
        self.pos = int(pos)
        self.ref = ref
        self.alts = alts.split(",")

        self.type = self.determine_variant_type()
        self.size = self.determine_variant_size()

        self.priorities = []

        for key, val in kwargs.items():
            if (not isinstance(key, str) or not isinstance(val, str)
                    or " " in key
                    or key.startswith("_")):
                continue
            self.__setattr__(key.lower(), val)

    def __str__(self):
        return f"{self.chrom}_{self.pos}_{self.ref}_{','.join(self.alts)}"

    def determine_variant_type(self):
        for allele in self.alts + [self.ref]:
            if len(allele) > 1:
                return "indel"
        return "snp"

    def determine_variant_size(self):
        return max(map(len, self.alts + [self.ref]))
# =============================================================================
#
#     def get_variant_chrom(self):
#         """Returns the chromosome name the variant is located on
#
#         Returns
#         -------
#         self.vcf_variant_chrom : str
#             Chromosome name the variant is located on
#         """
#         return self.vcf_variant_chrom
#
#     def set_variant_chrom(self, varchrom):
#         """Sets the chromosome name the variant is located on. Overwrites an already set chromosome name.
#
#         Parameters
#         ----------
#         varchrom : str
#             Chromosome name the variant is located on
#         :return:
#         """
#         self.vcf_variant_chrom = varchrom
#
#     def get_variant_pos(self):
#         """Returns the genomic position of the variant.
#
#         Returns
#         -------
#         self.vcf_variant_pos : nt
#             Genomic position of the variant
#         """
#         return self.vcf_variant_start
#
#     def set_variant_pos(self, varpos):
#         """Sets the variant genomic position. Overwrites an already set genomic position.
#
#         Parameters
#         ----------
#         varpos : int
#             Genomic position of the variant.
#         """
#         self.vcf_variant_start = varpos
#
#     def get_variant_ref_allele(self):
#         """Returns the reference allele of the variant.
#
#         Returns
#         -------
#         self.vcf_variant_ref : str
#             Reference allele of the variant
#         :return:
#         """
#         return self.vcf_variant_ref
#
#     def set_variant_ref_allele(self, varref):
#         """Sets the reference allele of the variant. Overwrites an already set reference allele.
#
#         Parameters
#         ----------
#         varref : str
#             Reference allele of the variant
#         """
#         self.vcf_variant_ref = varref
#
#     def get_variant_alt_alleles(self):
#         """Returns the alternative allele(s) of the variant.
#
#         Returns
#         -------
#         self.vcf_variant_alts : tuple
#             Alternative allele(s) of the variant
#         """
#         return self.vcf_variant_alts
#
#     def set_variant_alt_alleles(self, varalts):
#         """Sets the alternative allele(s) of the variants. Overwrites an already set alternative allele(s).
#
#         Parameters
#         ----------
#         varalts : str
#             Alternative allele(s) of the variant
#         """
#         self.vcf_variant_alts = tuple(varalts.split(","))
#
#     def get_variant_filter(self):
#         """Returns the filter (i.e. PASS) that was applied to the variant
#
#         Returns
#         -------
#         self.vcf_variant_filter : str
#
#         """
#         return self.vcf_variant_filter
#
#     def set_variant_filter(self, varfilter):
#         """Sets the filter applied to the variant. Overwrites any already set variant filter.
#
#         Parameters
#         ----------
#         varfilter : str
#
#         """
#         self.vcf_variant_filter = varfilter
#
#     def get_variant_type(self):
#         """Returns the variant type.
#
#         Returns
#         -------
#         self.vcf_variant_type : str
#             Variant type
#         """
#         return self.vcf_variant_type
#
#     def set_variant_type(self, vartype):
#         """Sets the type of the variant. Overwrites any already set variant type.
#
#         Parameters
#         ----------
#         vartype : str
#         """
#         self.vcf_variant_type = vartype
#
#     def get_variant_id(self):
#         """Constructs and returns a variant identifier.
#
#         The identifier is formed by combining the chromosome name and variant position, separated by an '_'
#         (i.e. CHROM_POS). This will be used as the identifier for contexts.
#
#         Returns
#         -------
#         str
#             Constructed variant identifier
#         """
#         return f"{self.vcf_variant_chrom}_{self.vcf_variant_start}"
#
#     def set_filter(self, filtername, filtervalue):
#         """Sets the value for a specified filter.
#
#         Parameters
#         ----------
#         filtername : str
#             Name of the used filter
#         filtervalue : str
#             Value for the used filter
#         """
#         self.priorityfilters[filtername] = filtervalue
#
#     def get_priority_filters(self):
#         """Returns the set priority filters.
#
#         Returns
#         -------
#         self.priorityfilters : dict
#             Map of set filters
#         """
#         return self.priorityfilters
#
#     def get_priority_filter(self, filtername):
#         """Returns the value of the specified priority filter.
#
#         Parameters
#         ----------
#         filtername : str
#             Name of the priority filter to get value of
#
#         Returns
#         -------
#         str or None
#             Value of the specified filter, None if the filter does not exist
#         """
#         if filtername in self.priorityfilters:
#             return self.priorityfilters[filtername]
#         return None
#
#     def set_priority_level(self, priority_filter, priority_level):
#         """Sets the priority level for a specified priority filter.
#
#         Parameters
#         ----------
#         priority_filter : str
#             Priority filter name
#         priority_level : int
#             Priority level to set for the variant
#         """
#         self.prioritylevels[priority_filter] = priority_level
#
#     def get_priority_level(self, filtername):
#         """Returns the determined priority level of the variant.
#
#         Parameters
#         ----------
#         filtername : str
#             Priority filter name
#
#         Returns
#         -------
#         str or None
#         """
#         if filtername in self.prioritylevels:
#             return self.prioritylevels[filtername]
#         return None
# =============================================================================

# =============================================================================
#     def to_string(self):
#         """Returns a String representation of the variant.
#
#         Returns
#         -------
#         str
#             String representation of the variant
#         """
#         return f"{self.vcf_variant_chrom}\t{self.vcf_variant_start}\t{self.vcf_variant_type}\t{self.vcf_variant_ref}" \
#             f"\t{self.vcf_variant_alts}\t{self.vcf_variant_filter}"
# =============================================================================
