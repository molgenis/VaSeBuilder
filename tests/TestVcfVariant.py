import unittest
from VcfVariant import VcfVariant


class TestVcfVariant(unittest.TestCase):
    def setUp(self):
        self.variant_chrom_answer = "1"
        self.variant_pos_answer = 100
        self.variant_ref_answer = "C"
        self.variant_alts_answer = ("A", "G")
        self.variant_filter_answer = ["PASS"]
        self.variant_type_answer = "snp"
        self.variant_id_answer = f"{self.variant_chrom_answer}_{self.variant_pos_answer}"
        self.vcfvariant_obj_answer = VcfVariant(self.variant_chrom_answer, self.variant_pos_answer,
                                                self.variant_ref_answer, self.variant_alts_answer,
                                                self.variant_filter_answer, self.variant_type_answer)
        self.to_string_answer = f"{self.variant_chrom_answer}\t{self.variant_pos_answer}\t{self.variant_type_answer}" \
                                f"\t{self.variant_ref_answer}\t{self.variant_alts_answer}\t{self.variant_filter_answer}"

    def test_get_variant_chrom(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_chrom(), self.variant_chrom_answer, "Both variant "
                         f"chromosomes should have been {self.variant_chrom_answer}")

    def test_set_variant_chrom(self):
        varchrom_toset = "21"
        self.vcfvariant_obj_answer.set_variant_chrom(varchrom_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_chrom, varchrom_toset, "Both set variant chromosomes "
                         f"should have been {varchrom_toset}")

    def test_get_variant_pos(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_pos(), self.variant_pos_answer, "Both variant positions"
                         f" should have been {self.variant_pos_answer}")

    def test_set_variant_pos(self):
        varpos_toset = 25000
        self.vcfvariant_obj_answer.set_variant_pos(varpos_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_start, varpos_toset, f"Both set variant positions "
                         f"should have been {varpos_toset}")

    def test_get_variant_ref_allele(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_ref_allele(), self.variant_ref_answer, "Both variant "
                         f"reference alleles should have been: {self.variant_ref_answer}")

    def test_set_variant_ref_allele(self):
        varref_toset = "A"
        self.vcfvariant_obj_answer.set_variant_ref_allele(varref_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_ref, varref_toset, "The set variant reference allele "
                         f"should have been: {varref_toset}")

    def test_get_variant_alt_alleles(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_alt_alleles(), self.variant_alts_answer, "Both variant "
                         f"alternative alleles should have been: {self.variant_alts_answer}")

    def test_set_variant_alt_alleles(self):
        varalts_toset = ("T", "G")
        self.vcfvariant_obj_answer.set_variant_alt_alleles(varalts_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_alts, varalts_toset, "The set variant alternative "
                         f"alleles should have been: {varalts_toset}")

    def test_get_variant_filter(self):
        self.assertListEqual(self.vcfvariant_obj_answer.get_variant_filter(), self.variant_filter_answer, "Both variant"
                             f" filters should have been: {self.variant_filter_answer}")

    def test_set_variant_filter(self):
        varfilter_touse = ['q10']
        self.vcfvariant_obj_answer.set_variant_filter(varfilter_touse)
        self.assertListEqual(self.vcfvariant_obj_answer.get_variant_filter(), varfilter_touse, "The set variant filter "
                             f"should have been: {varfilter_touse}")

    def test_get_variant_type(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_type(), self.variant_type_answer, "Both variant types "
                         f"should have been: {self.variant_type_answer}")

    def test_set_variant_type(self):
        vartype_touse = "indel"
        self.vcfvariant_obj_answer.set_variant_type(vartype_touse)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_type, vartype_touse, "The set variant type should have "
                         f"been {vartype_touse}")

    def test_get_variant_id(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_id(), self.variant_id_answer, "Both variant identifiers"
                         f" should have been: {self.variant_id_answer}")

    def test_to_string(self):
        self.assertEqual(self.vcfvariant_obj_answer.to_string(), self.to_string_answer, f"Both VcfVariant string "
                         f"representations should have been: {self.to_string_answer}")
