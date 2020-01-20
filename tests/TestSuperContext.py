import unittest
from SuperContext import SuperContext
from VariantContext import VariantContext


class TestSuperContext(unittest.TestCase):
    def setUp(self):
        # Create some stuff
        self.chrom_answer = "21"
        self.origin_answer = 150
        self.start_answer = 100
        self.end_answer = 400
        self.context_answer = []
        self.genomic_region = ""

        self.first_varcon = VariantContext("21_150", "SAMPLE_A", self.chrom_answer, self.origin_answer,
                                           self.start_answer, 250, [], [])
        self.second_varcon = VariantContext("21_300", "SAMPLE_A", self.chrom_answer, 300, 200, self.end_answer, [], [])
        self.super_context = SuperContext(self.first_varcon, self.second_varcon)

    # Tests that the correct chromosome is returned (and therefore saved)
    def test_get_chrom(self):
        self.assertEqual(self.super_context.get_chrom(), self.chrom_answer,
                         f"The super context chromosome should have been {self.chrom_answer}")

    # Tests that the origin position is returned (and therefore saved)
    def test_get_origin(self):
        self.assertEqual(self.super_context.get_origin(), self.origin_answer,
                         f"the super context origin position should have been {self.origin_answer}")

    # Tests that the correct start position is returned (and therefore saved)
    def test_get_start(self):
        self.assertEqual(self.super_context.get_start(), self.start_answer,
                         f"The super context start position should have been {self.start_answer}")

    # Tests that the correct end position is returned (and therefore saved)
    def test_get_end(self):
        self.assertEqual(self.super_context.get_end(), self.end_answer,
                         f"The super context end position should have been {self.end_answer}")

    # Tests that the correct variant contexts are returned.
    def test_get_variant_contexts(self):
        varcons_answer = [self.first_varcon.get_context(), self.second_varcon.get_context()]
        varcons_received = [x.get_context() for x in self.super_context.get_variant_contexts()]
        self.assertListEqual(varcons_received, varcons_answer,
                             f"The received variant contexts are not what should have been.")

    # Tests that the super context has the correct number of variant contexts
    def test_get_num_of_variant_contexts(self):
        self.assertEqual(self.super_context.get_num_of_variant_contexts(), 2,
                         f"There should have been 2 variant contexts in the super context")

    def test_get_context(self):
        self.assertListEqual(self.super_context.get_context(), self.context_answer,
                             f"The returned super context window should have been {self.context_answer}")

    def test_get_genomic_region(self):
        self.assertEqual(self.super_context.get_genomic_region(), self.genomic_region,
                         f"The returned genomic region should have been {self.genomic_region}")

    # Tests thst a new variant context is added.
    def test_add_variant_context(self):
        varcon_to_add = VariantContext("21_450", "SAMPLE_A", "21", 450, 400, 500, [], [])
        context_answer = varcon_to_add.get_context()

        self.super_context.add_variant_context(varcon_to_add)
        added_context = self.super_context.get_variant_contexts()[-1].get_context()
        self.assertListEqual(added_context, context_answer,
                             f"The added context should have been {context_answer}")

    def test_add_variant_contexts(self):

    # def test_add_variant_contexts(self):

    # Tests that the a new, smaller start position is set as the super context start position
    def test_determine_super_start_smallerstart(self):
        smaller_start = 75
        self.super_context.determine_super_start(smaller_start)
        self.assertEqual(self.super_context.get_start(), smaller_start,
                         f"The new super context start position should have been {smaller_start}")

    # Tests that the new, larger start position is not set as the super context start position
    def test_determine_super_start_largerstart(self):
        larger_start = 125
        self.super_context.determine_super_start(larger_start)
        self.assertNotEqual(self.super_context.get_start(), larger_start,
                            f"The super context start position should not have been {larger_start}")

    # Tests that the new, larger end position is set as the new super context end position
    def test_determine_super_end_largerend(self):
        larger_end = 450
        self.super_context.determine_super_end(larger_end)
        self.assertEqual(self.super_context.get_end(), larger_end,
                         f"The super context end position should have been {larger_end}")

    # Tests that the enw, smaller end position is not set as the new super context end position
    def test_determine_super_end_smallerend(self):
        smaller_end = 350
        self.super_context.determine_super_end(smaller_end)
        self.assertNotEqual(self.super_context.get_end(), smaller_end,
                            f"The super context end position should not have been {smaller_end}")

    # Tests that the default super context indeed has no gaps.
    def test_has_gaps_nogaps(self):
        self.assertFalse(self.super_context.has_gaps(), "The default super context should not have any gaps")

    # Tests that the super context has a gap after adding a non overlapping variant context
    def test_has_gaps_onegap(self):
        non_overlap_varcon = VariantContext("21_1000", "SAMPLE_A", "21", 900, 1000, 1250, [], [])
        self.super_context.add_variant_context(non_overlap_varcon)
        self.assertTrue(self.super_context.has_gaps(),
                        "The super context with the added variant context should have had a gap.")

    # def test_determine_rightmost_context(self):

    def test_varcons_overlap_positive(self):
        self.assertTrue(self.super_context.varcons_overlap(self.first_varcon, self.second_varcon),
                        "The two variant context should overlap.")

    def test_varcons_overlap_negative(self):
        non_overlap_varcon = VariantContext("21_1000", "SAMPLE_A", "21", 900, 1000, 1250, [], [])
        self.assertFalse(self.super_context.varcons_overlap(self.first_varcon, non_overlap_varcon),
                         "The two variant contexts should not overlap")

    # def determine_largest_variantcontext(self):
    # def test_get_varcons_by_leftpos(self):
    # def fix_gaps(self):
