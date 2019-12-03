import unittest
from VcfVariant import VcfVariant


class TestVcfVariant(unittest.TestCase):
    def setUp(self):
        self.variant_chrom_answer = "1"
        self.variant_pos_answer = 100
        self.variant_ref_answer = "C"
        self.variant_alts_answer = ("A", "G")
        self.variant_filter_answer = "PASS"
        self.variant_type_answer = "snp"
        self.variant_id_answer = f"{self.variant_chrom_answer}_{self.variant_pos_answer}"
        self.priority_filter_col = "Relevance"
        self.priority_filter_value = "BENIGN"
        self.priority_filter_answer = {self.priority_filter_col: self.priority_filter_value}
        self.priority_level_answer = 2
        self.vcfvariant_obj_answer = VcfVariant(self.variant_chrom_answer, self.variant_pos_answer,
                                                self.variant_ref_answer, self.variant_alts_answer)
        self.to_string_answer = f"{self.variant_chrom_answer}\t{self.variant_pos_answer}\t{self.variant_type_answer}" \
                                f"\t{self.variant_ref_answer}\t{self.variant_alts_answer}\t{self.variant_filter_answer}"

    # Tests that the correct variant chromosome name is returned
    def test_get_variant_chrom(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_chrom(), self.variant_chrom_answer, "Both variant "
                         f"chromosomes should have been {self.variant_chrom_answer}")

    # Tests that the variant chromosome is indeed set correctly.
    def test_set_variant_chrom(self):
        varchrom_toset = "21"
        self.vcfvariant_obj_answer.set_variant_chrom(varchrom_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_chrom, varchrom_toset, "Both set variant chromosomes "
                         f"should have been {varchrom_toset}")

    # Tests that the correct variant position is returned.
    def test_get_variant_pos(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_pos(), self.variant_pos_answer, "Both variant positions"
                         f" should have been {self.variant_pos_answer}")

    # Tests that the variant position is indeed set.
    def test_set_variant_pos(self):
        varpos_toset = 25000
        self.vcfvariant_obj_answer.set_variant_pos(varpos_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_start, varpos_toset, f"Both set variant positions "
                         f"should have been {varpos_toset}")

    # Tests that the correct variant reference allele is returned.
    def test_get_variant_ref_allele(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_ref_allele(), self.variant_ref_answer, "Both variant "
                         f"reference alleles should have been: {self.variant_ref_answer}")

    # Tests that the variant reference allele is indeed set.
    def test_set_variant_ref_allele(self):
        varref_toset = "A"
        self.vcfvariant_obj_answer.set_variant_ref_allele(varref_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_ref, varref_toset, "The set variant reference allele "
                         f"should have been: {varref_toset}")

    # Tests that the correct alternative alleles are returned.
    def test_get_variant_alt_alleles(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_alt_alleles(), self.variant_alts_answer, "Both variant "
                         f"alternative alleles should have been: {self.variant_alts_answer}")

    # Tests that the string variant alternative alleles are correctly set as a tuple
    def test_set_variant_alt_alleles(self):
        varalts_answer = ("T", "G")
        varalts_toset = "T, G"
        self.vcfvariant_obj_answer.set_variant_alt_alleles(varalts_toset)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_alts, varalts_answer, "The set variant alternative "
                         f"alleles should have been: {varalts_answer}")

    # Tests that the correct variant filter is returned.
    def test_get_variant_filter(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_filter(), self.variant_filter_answer,
                         f"Both variant filters should have been: {self.variant_filter_answer}")

    # Tests that the variant filter is indeed set.
    def test_set_variant_filter(self):
        varfilter_touse = "q10"
        self.vcfvariant_obj_answer.set_variant_filter(varfilter_touse)
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_filter(), varfilter_touse, "The set variant filter "
                         f"should have been: {varfilter_touse}")

    # Tests that the correct variant type is returned.
    def test_get_variant_type(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_type(), self.variant_type_answer, "Both variant types "
                         f"should have been: {self.variant_type_answer}")

    # Tests that the variant type is indeed set.
    def test_set_variant_type(self):
        vartype_touse = "indel"
        self.vcfvariant_obj_answer.set_variant_type(vartype_touse)
        self.assertEqual(self.vcfvariant_obj_answer.vcf_variant_type, vartype_touse, "The set variant type should have "
                         f"been {vartype_touse}")

    # Tests that the variant ID is returned as CHROM_POS
    def test_get_variant_id(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_variant_id(), self.variant_id_answer, "Both variant identifiers"
                         f" should have been: {self.variant_id_answer}")

    # Tests that a filter is indeed set correctly.
    def test_set_filter(self):
        priority_filter_name = "noot"
        priority_filter_answer = "mies"
        self.vcfvariant_obj_answer.set_filter("noot", "mies")
        self.assertEqual(self.vcfvariant_obj_answer.priorityfilters[priority_filter_name], priority_filter_answer,
                         f"The set priority filter should have been {priority_filter_answer}")

    # Tests that the correct priority filters are returned.
    def test_get_priority_filters(self):
        self.assertDictEqual(self.vcfvariant_obj_answer.get_priority_filters(), self.priority_filter_answer,
                             f"The set priority filters should have been {self.priority_filter_answer}")

    # Tests that the correct priority filter is returned.
    def test_get_priority_filter(self):
        self.assertEqual(self.vcfvariant_obj_answer.get_priority_filter(self.priority_filter_col),
                         self.priority_filter_value,
                         f"The value for {self.priority_filter_col} should have been {self.priority_filter_value}")

    # Test that None is returned when requesting the value of an incorrect filter name.
    def test_get_priority_filter_incorrectfilter(self):
        self.assertIsNone(self.vcfvariant_obj_answer.get_priority_filter("aap"),
                          "Filter aap should not have been a valid filter and return None.")

    # Tests that the expected string representation is returned.
    def test_to_string(self):
        self.assertEqual(self.vcfvariant_obj_answer.to_string(), self.to_string_answer, f"Both VcfVariant string "
                         f"representations should have been: {self.to_string_answer}")
