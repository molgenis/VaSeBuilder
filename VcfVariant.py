"""Object to save the essential data of a VCF variant"""


class VcfVariant:
    """Saves necessary data of a genomic variant from a VCF/BCF file.

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
    def __init__(self, varchrom, varstart, varref, varalts, varfilter, vartype):
        self.vcf_variant_chrom = varchrom
        self.vcf_variant_start = varstart
        self.vcf_variant_ref = varref
        self.vcf_variant_alts = varalts
        self.vcf_variant_filter = varfilter
        self.vcf_variant_type = vartype

    def get_variant_chrom(self):
        """Returns the chromosome name the variant is located on

        Returns
        -------
        self.vcf_variant_chrom : str
            Chromosome name the variant is located on
        """
        return self.vcf_variant_chrom

    def set_variant_chrom(self, varchrom):
        """Sets the chromosome name the variant is located on. Overwrites an already set chromosome name.

        Parameters
        ----------
        varchrom : str
            Chromosome name the variant is located on
        :return:
        """
        self.vcf_variant_chrom = varchrom

    def get_variant_pos(self):
        """Returns the genomic position of the variant.

        Returns
        -------
        self.vcf_variant_pos : nt
            Genomic position of the variant
        """
        return self.vcf_variant_start

    def set_variant_pos(self, varpos):
        """Sets the variant genomic position. Overwrites an already set genomic position.

        Parameters
        ----------
        varpos : int
            Genomic position of the variant.
        """
        self.vcf_variant_start = varpos

    def get_variant_ref_allele(self):
        """Returns the reference allele of the variant.

        Returns
        -------
        self.vcf_variant_ref : str
            Reference allele of the variant
        :return:
        """
        return self.vcf_variant_ref

    def set_variant_ref_allele(self, varref):
        """Sets the reference allele of the variant. Overwrites an already set reference allele.

        Parameters
        ----------
        varref : str
            Reference allele of the variant
        """
        self.vcf_variant_ref = varref

    def get_variant_alt_alleles(self):
        """Returns the alternative allele(s) of the variant.

        Returns
        -------
        self.vcf_variant_alts : tuple
            Alternative allele(s) of the variant
        """
        return self.vcf_variant_alts

    def set_variant_alt_alleles(self, varalts):
        """Sets the alternative allele(s) of the variants. Overwrites an already set alternative allele(s).

        Parameters
        ----------
        varalts : tuple
            Alternative allele(s) of the variant
        """
        self.vcf_variant_alts = varalts

    def get_variant_filter(self):
        """Returns the filter (i.e. PASS) that was applied to the variant

        Returns
        -------
        self.vcf_variant_filter : str

        """
        return self.vcf_variant_filter

    def set_variant_filter(self, varfilter):
        """Sets the filter applied to the variant. Overwrites any already set variant filter.

        Parameters
        ----------
        varfilter : str

        """
        self.vcf_variant_filter = varfilter

    def get_variant_type(self):
        """Returns the variant type.

        Returns
        -------
        self.vcf_variant_type : str
            Variant type
        """
        return self.vcf_variant_type

    def set_variant_type(self, vartype):
        """Sets the type of the variant. Overwrites any already set variant type.

        Parameters
        ----------
        vartype : str
        """
        self.vcf_variant_type = vartype

    def get_variant_id(self):
        """Constructs and returns a variant identifier.

        The identifier is formed by combining the chromosome name and variant position, separated by an '_'
        (i.e. CHROM_POS). This will be used as the identifier for contexts.

        Returns
        -------
        str
            Constructed variant identifier
        """
        return f"{self.vcf_variant_chrom}_{self.vcf_variant_start}"

    def to_string(self):
        """Returns a String representation of the variant.

        Returns
        -------
        str
            String representation of the variant
        """
        return f"{self.vcf_variant_chrom}\t{self.vcf_variant_start}\t{self.vcf_variant_type}\t{self.vcf_variant_ref}" \
            f"\t{self.vcf_variant_alts}\t{self.vcf_variant_filter}"
