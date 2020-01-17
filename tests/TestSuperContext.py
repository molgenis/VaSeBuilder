import unittest
from SuperContext import SuperContext
from VariantContext import VariantContext


class TestSuperContext(unittest.TestCase):
    def setUp(self):
        # Create some stuff
        self.chrom_answer = ""
        self.origin_answer = 0
        self.start_answer = 0
        self.end_answer = 0
        self.context_answer = []
        self.genomic_region = ""

        self.first_varcon = VariantContext()
        self.second_varcon = VariantContext()
        self.super_context = SuperContext(self.first_varcon, self.second_varcon)

    def test_get_chrom(self):
        self.assertEqual(self.super_context.get_chrom(), self.chrom_answer,
                         f"The super context chromosome should have been {self.chrom_answer}")

    def test_get_origin(self):
        self.assertEqual(self.super_context.get_origin(), self.origin_answer,
                         f"the super context origin position should have been {self.origin_answer}")

    def test_get_start(self):
        self.assertEqual(self.super_context.get_start(), self.start_answer,
                         f"The super context start position should have been {self.start_answer}")

    def test_get_end(self):
        self.assertEqual(self.super_context.get_end(), self.end_answer,
                         f"The super context end position should have been {self.end_answer}")

    def test_get_variant_contexts(self):
        varcons_answer = [self.first_varcon.get_context(), self.second_varcon.get_context()]
        varcons_received = [x.get_context() for x in self.super_context.get_variant_contexts()]
        self.assertListEqual(varcons_received, varcons_answer,
                             f"The received variant contexts are not what should have been.")

    def test_get_num_of_variant_contexts(self):
        self.assertEqual(self.super_context.get_num_of_variant_contexts(), 2,
                         f"")

    def test_get_context(self):
        self.assertListEqual(self.super_context.get_context(), self.context_answer,
                             f"The returned super context window should have been {self.context_answer}")

    def test_get_genomic_region(self):
        self.assertEqual(self.super_context.get_genomic_region(), self.genomic_region,
                         f"The returned genomic region should have been {self.genomic_region}")

    # def test_add_variant_context(self):
    # def test_add_variant_contexts(self)

