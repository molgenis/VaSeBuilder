import unittest
from CompoundContext import CompoundContext
from VariantContext import VariantContext


class TestCompoundContext(unittest.TestCase):
    def setUp(self):
        self.chr_a_answer = "21"
        self.chr_b_answer = "B"
        self.sample_answer = "SAMPLE_AB"
        self.a_start_answer = 100
        self.b_start_answer = 250
        self.a_origin_answer = 150
        self.b_origin_answer = 300
        self.a_end_answer = 200
        self.b_end_answer = 350

        self.variant_context_a = VariantContext(f"{self.chr_a_answer}_{self.a_origin_answer}", self.sample_answer,
                                                self.chr_a_answer, self.a_origin_answer, self.a_start_answer,
                                                self.a_end_answer, [], [])
        self.variant_context_b = VariantContext(f"{self.chr_b_answer}_{self.b_origin_answer}", self.sample_answer,
                                                self.chr_b_answer, self.b_origin_answer, self.b_start_answer,
                                                self.b_end_answer, [], [])
        self.compound_context = CompoundContext(self.variant_context_a, self.variant_context_b)

    # Tests that the correct variant contexts are returned
    def test_get_variants_contexts(self):
        self.contexts_answer = []

    # Tests that the the CompoundContext has the correct number of VariantContext
    def test_get_num_of_variant_contexts(self):
        self.assertEqual(self.compound_context.get_num_of_variant_contexts(), 2,
                         "The number of variant contexts in the compound context should have been 2")

    # Tests that a None variant context will indeed not be added
    def test_add_variant_context_none(self):
        self.compound_context.add_variant_context(None)
        self.assertEqual(self.compound_context.get_num_of_variant_contexts(), 2,
                         "The number of variant contexts still be 2 as None variant contexts will not be added")

    # Tests that a variant context can indeed be added.
    def test_add_variant_context(self):
        varcon_to_add = VariantContext("6_550", "SAMPLE_A", "6", 550, 500, 600, [], [])
        context_answer = ["6", 550, 500, 600]
        self.compound_context.add_variant_context(varcon_to_add)
        received_context = self.compound_context.get_variant_contexts()[-1]
        self.assertListEqual(received_context, context_answer,
                             f"The context {context_answer} should have been added")

    # Tests that a None is not added to the list of variant contexts
    def test_add_variant_contexts_none(self):
        self.compound_context.add_variant_contexts(None)
        self.assertEqual(self.compound_context.get_num_of_variant_contexts(), 2,
                         "The number of variant contexts should have remained 2")

    # Tests that a list of Nones is not added to the compound context
    def test_add_variant_contexts_nonelist(self):
        varcons_to_add = [None, None]
        self.compound_context.add_variant_contexts(varcons_to_add)
        self.assertEqual(self.compound_context.get_num_of_variant_contexts(), 2,
                         "The number of variant contexts should have remained 2")

    # Tests that multiple variants contexts can be added
    def test_add_variant_contexts(self):
        varcon_to_add_1 = VariantContext("1_150", "SAMPLE_A", "1", 150, 100, 200, [], [])
        varcon_to_add_2 = VariantContext("2_250", "SAMPLE_A", "2", 250, 200, 300, [], [])
        varcons_toadd = [varcon_to_add_1, varcon_to_add_2]
        self.compound_context.add_variant_contexts(varcons_toadd)
        self.assertEqual(self.compound_context.get_num_of_variant_contexts(), 4,
                         "The compound context should have 4 variant contexts as we wanted to add two extra")

    # def test_determine_priority_level(self):
    # def test_remove_variant_context(self):
