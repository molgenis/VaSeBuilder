"""Object to save the essential data of a VCF variant"""


class VcfVariant:
    def __init__(self, varchrom, varstart, varref, varalts, varfilter, vartype):
        self.vcf_variant_chrom = varchrom
        self.vcf_variant_start = varstart
        self.vcf_variant_ref = varref
        self.vcf_variant_alts = varalts
        self.vcf_variant_filter = varfilter
        self.vcf_variant_type = vartype

    # Returns the chromosome name of the variant
    def get_variant_chrom(self):
        return self.vcf_variant_chrom

    # Returns the VCF variant start position
    def get_variant_pos(self):
        return self.vcf_variant_start

    # Returns the VCF variant reference allele
    def get_variant_ref_allele(self):
        return self.vcf_variant_ref

    # Returns the VCF variant alternative alleles
    def get_variant_alt_alleles(self):
        return self.vcf_variant_alts

    # Returns the VCF variant filter
    def get_variant_filter(self):
        return self.vcf_variant_filter

    # Returns the VCF variant type (SNP/Indel, etc)
    def get_variant_type(self):
        return self.vcf_variant_type

    # Returns a variant identifier in the form CHROM_POS
    def get_variant_id(self):
        return f"{self.vcf_variant_chrom}_{self.vcf_variant_start}"

    # ToString method
    def to_string(self):
        return f"{self.vcf_variant_chrom}\t{self.vcf_variant_start}\t{self.vcf_variant_type}\t{self.vcf_variant_ref}" \
            f"\t{self.vcf_variant_alts}\t{self.vcf_variant_filter}"
