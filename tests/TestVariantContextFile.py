import unittest

# Import the required VaSe classes
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile


class TestVariantContextFile(unittest.TestCase):
    def setUp(self):
        # Construct the bam reads to use
        self.context_chrom_answer = '21'

        self.va_read_id = "vaRead1"
        self.vd_read_id = "vdRead1"
        self.a_read_id = "aRead1"
        self.d_read_id = "dRead1"
        self.a_read_start_pos_1 = 9411000
        self.a_read_start_pos_2 = 9411200
        self.a_read_end_pos = 9411350
        self.d_read_start_pos_1 = 9411150
        self.d_read_start_pos_2 = 9411350
        self.d_read_end_pos = 9411500
        self.bamread_len = 151
        self.a_read1_seq = "TTTAGATGGG"
        self.a_read2_seq = "ATTTCTAGTT"
        self.d_read1_seq = "AGAAAAAGTC"
        self.d_read2_seq = "TGCCTTTTCA"
        self.a_read1_quals = "=====<===="
        self.a_read2_quals = ">???>?????"
        self.d_read1_quals = "><=???>==<"
        self.d_read2_quals = "TGCCTTTTCA"
        self.read_mapq = 40
        self.a_read_1 = DonorBamRead(self.a_read_id, "1", self.context_chrom_answer, self.a_read_start_pos_1,
                                     self.bamread_len, self.a_read1_seq, self.a_read1_quals, self.read_mapq)
        self.a_read_2 = DonorBamRead(self.a_read_id, "2", self.context_chrom_answer, self.a_read_start_pos_2,
                                     self.bamread_len, self.a_read2_seq, self.a_read2_quals, self.read_mapq)
        self.d_read_1 = DonorBamRead(self.d_read_id, "1", self.context_chrom_answer, self.d_read_start_pos_1,
                                     self.bamread_len, self.d_read1_seq, self.d_read1_quals, self.read_mapq)
        self.d_read_2 = DonorBamRead(self.d_read_id, "2", self.context_chrom_answer, self.d_read_start_pos_2,
                                     self.bamread_len, self.d_read2_seq, self.d_read2_quals, self.read_mapq)
        self.va_read_1 = DonorBamRead(self.va_read_id, "1", self.context_chrom_answer, self.a_read_start_pos_1,
                                      self.bamread_len, self.a_read1_seq, self.a_read1_quals, self.read_mapq)
        self.va_read_2 = DonorBamRead(self.va_read_id, "2", self.context_chrom_answer, self.a_read_start_pos_2,
                                      self.bamread_len, self.a_read2_seq, self.a_read2_quals, self.read_mapq)
        self.vd_read_1 = DonorBamRead(self.vd_read_id, "1", self.context_chrom_answer, self.d_read_start_pos_1,
                                      self.bamread_len, self.d_read1_seq, self.d_read1_quals, self.read_mapq)
        self.vd_read_2 = DonorBamRead(self.vd_read_id, "2", self.context_chrom_answer, self.d_read_start_pos_2,
                                      self.bamread_len, self.d_read1_seq, self.d_read1_quals, self.read_mapq)
        
        # Create the variables containing context answers and values
        self.context_id_answer = "21_9411259"
        self.context_sample_answer = "testsample"
        self.context_origin_answer = 9411250
        self.acc_context_start_answer = 9411000
        self.don_context_start_answer = 9411150
        self.var_context_start_answer = self.acc_context_start_answer
        self.acc_context_end_answer = 9411350
        self.don_context_end_answer = 9411500
        self.var_context_end_answer = self.don_context_end_answer
        self.acc_context_reads_answer = [self.a_read_1, self.a_read_2]
        self.don_context_reads_answer = [self.d_read_1, self.d_read_2]
        self.var_context_a_reads_answer = [self.va_read_1, self.va_read_2]
        self.var_context_d_reads_answer = [self.vd_read_1, self.vd_read_2]
        self.acceptor_context_answer = OverlapContext(self.context_id_answer, self.context_sample_answer,
                                                      self.context_chrom_answer, self.context_origin_answer,
                                                      self.acc_context_start_answer, self.acc_context_end_answer,
                                                      self.acc_context_reads_answer)
        self.donor_context_answer = OverlapContext(self.context_id_answer, self.context_sample_answer,
                                                   self.context_chrom_answer, self.context_origin_answer,
                                                   self.don_context_start_answer, self.don_context_end_answer,
                                                   self.don_context_reads_answer)
        
        # Create the variables containing all the VariantContext answers
        self.variantContextAnswer = VariantContext(self.context_id_answer, self.context_sample_answer,
                                                   self.context_chrom_answer, self.context_origin_answer,
                                                   self.var_context_start_answer, self.var_context_end_answer,
                                                   self.var_context_a_reads_answer, self.var_context_d_reads_answer)

        self.setAccContextAnswer = OverlapContext('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050,
                                                  self.acc_context_reads_answer)
        self.setDonContextAnswer = OverlapContext('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100,
                                                  self.don_context_reads_answer)
        self.setVarContextAnswer = VariantContext('21_9411000', 'testanswer', '21', 9411000, 9410900, 9411100,
                                                  self.var_context_a_reads_answer, self.var_context_d_reads_answer,
                                                  self.setAccContextAnswer, self.setDonContextAnswer)

        self.varcon_fields_answer = {1: "variant context id", 2: "sample id", 3: "chromosome", 4: "origin",
                                     5: "start pos", 6: "end pos", 7: "acceptor context", 8: "donor context",
                                     9: "number of acceptor reads", 10: "number of donor reads",
                                     11: "acceptor/donor ratio", 12: "acceptor read ids", 13: "donor read ids"}

        # Create the variables containg the VariantContextFile answers
        self.filterListToUse = ["aap", "noot", "mies"]
        
        # Create the variables containing the answer of the VariantContextFile
        self.varcon_file = VariantContextFile()
        self.varcon_file.set_variant_context(self.context_id_answer, self.variantContextAnswer)
        self.varcon_file.set_acceptor_context(self.context_id_answer, self.acceptor_context_answer)
        self.varcon_file.set_donor_context(self.context_id_answer, self.donor_context_answer)
        
        # Construct the other VariantContextFile objects to use for the set operation tests
        self.pos_otherVariantContextFile = VariantContextFile()
        self.pos_otherVariantContextFile.set_variant_context(self.context_id_answer, self.variantContextAnswer)
        self.neg_otherVariantContextFile = VariantContextFile()
        self.neg_otherVariantContextFile.add_variant_context('20_150', 'testsample2', '20', 150, 0, 350,
                                                             self.var_context_d_reads_answer,
                                                             self.var_context_a_reads_answer)

        self.read_varcontxt_loc = "testdata/varcon.txt"
        self.read_varcontxt_entry1 = "21_9411327\tSAM001\t21\t9411327\t9411192\t9411908\t716\t716\t3\t3\t1.0\t" \
                                     "aRead11;aRead12;aRead13\tdRead11;dRead12;dRead13"
        self.read_varcontxt_entry2 = "22_9600250\tSAM002\t22\t9600250\t9600100\t9600500\t400\t400\t3\t3\t1.0\t" \
                                     "aRead21;aRead22;aRead23\tdRead21;dRead22;dRead23"
        self.read_varcontxt_entry3 = "21_9900000\tSAM003\t21\t9900000\t9899900\t9900250\t350\t350\t3\t3\t1.0\t" \
                                     "aRead31;aRead32;aRead33\tdRead31;dRead32;dRead33"
        self.acccontxt_loc = "testdata/acceptorcontexts.txt"
        self.doncontxt_loc = "testdata/donorcontexts.txt"

    # ====================TESTS FOR THE VARIANT CONTEXT WITHIN A VARIANT CONTEXT FILE====================
    # Tests that all the correct contexts are returned
    def test_get_variant_contexts(self):
        variant_contexts_answer = [self.variantContextAnswer.to_string()]
        obtained_contexts_answer = [x.to_string() for x in self.varcon_file.get_variant_contexts()]
        self.assertListEqual(obtained_contexts_answer, variant_contexts_answer,
                             "The returned variant contexts are not what was expected")

    # Tests that the variant context file has an existing variant context
    def test_has_variant_context_pos(self):
        self.assertTrue(self.varcon_file.has_variant_context(self.context_id_answer), "Variant context "
                        f"{self.context_id_answer} should have been in the variant context file.")

    # Tests that the variant context file does not have a non existing variant context
    def test_has_variant_context_neg(self):
        nonexisting_context_id = "1_100"
        self.assertFalse(self.varcon_file.has_variant_context(nonexisting_context_id), "There should not have "
                         f"been a variant context with identifier {nonexisting_context_id}")

    # Tests that an existing variant context is returned.
    def test_get_variant_context(self):
        self.assertEqual(self.varcon_file.get_variant_context(self.context_id_answer).to_string(),
                         self.variantContextAnswer.to_string(), "The returned variant context is not what was expected")

    # Tests that None is returned as a non existing variant context is specified
    def test_get_variant_context_none(self):
        self.assertIsNone(self.varcon_file.get_variant_context('22_9411255'),
                          "The requested variant context should not have existed and should have therefore been None")

    # Tests that the correct number of variant contexts
    def test_get_number_of_contexts(self):
        self.assertEqual(self.varcon_file.get_number_of_contexts(), 1, "There should have been only 1 variant "
                                                                              "context")

    # Tests that the correct list of variant context identifiers are returned
    def get_variant_context_ids(self):
        varcon_ids_answer = []
        self.assertListEqual(self.varcon_file.get_variant_context_ids(), varcon_ids_answer,
                             f"The returned list of varian context identifiers should have been {varcon_ids_answer}")

    # Tests that the correct variant context field map is returned
    def test_get_variant_context_fields(self):
        self.assertDictEqual(self.varcon_file.get_variant_context_fields(), self.varcon_fields_answer,
                             f"The returned variant context field should have been {self.varcon_fields_answer}")

    # Tests that all the acceptor reads for an existing variant context are returned
    def test_get_variant_context_areads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_variant_context_areads(self.context_id_answer)]
        answer_reads = [x.to_string() for x in self.var_context_a_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, f"The returned reads should have been {answer_reads}")

    # Tests that the correct variant context acceptor reads are returned
    def test_get_all_variant_context_acceptor_reads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_all_variant_context_acceptor_reads()]
        answer_reads = [x.to_string() for x in self.var_context_a_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same "
                                                           "data")

    # Tests that the correct variant context donor reads are returned
    def test_get_variant_context_dreads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_variant_context_dreads(self.context_id_answer)]
        answer_reads = [x.to_string() for x in self.var_context_d_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, f"The returned reads should have been {answer_reads}")

    # Tests that the correct variant context donor reads are returned
    def test_get_all_variant_context_donor_reads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_all_variant_context_donor_reads()]
        answer_reads = [x.to_string() for x in self.var_context_d_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same "
                                                           "data")

    # Tests that the correct variant context acceptor read identifiers are returned
    def test_get_all_variant_context_acceptor_read_ids(self):
        vca_read_ids = [self.va_read_id]
        self.assertListEqual(self.varcon_file.get_all_variant_context_acceptor_read_ids(), vca_read_ids,
                             f"The list of returned acceptor read ids should have been: {vca_read_ids}")

    # Tests that the correct variant context donor read identifiers are returned
    def test_get_all_variant_context_donor_read_ids(self):
        vcd_read_ids = [self.vd_read_id]
        self.assertListEqual(self.varcon_file.get_all_variant_context_donor_read_ids(), vcd_read_ids,
                             f"The list of returned donor read ids should have been: {vcd_read_ids}")

    # ===============TESTS FOR THE ACCEPTOR CONTEXT WITHIN A VARIANT CONTEXT FILE===============
    # Tests that the correct acceptor context is returned
    def test_get_acceptor_context(self):
        self.assertEqual(self.varcon_file.get_acceptor_context(self.context_id_answer).to_string(),
                         self.acceptor_context_answer.to_string(),
                         "The returned acceptor contexts is not what was expected")

    # Tests that None is returned as a non existing acceptor context is specified
    def test_get_acceptor_context_none(self):
        self.assertIsNone(self.varcon_file.get_acceptor_context("22_9411255"),
                          "The requested acceptor context should not have existed and should have therefore been None")

    # Tests that the correct acceptor reads are returned
    def test_get_acceptor_context_reads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_acceptor_context_reads(self.context_id_answer)]
        answer_reads = [x.to_string() for x in self.acc_context_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, "The returned acceptor context reads should have been "
                             f"{answer_reads}")

    # ===============TESTS FOR THE DONOR CONTEXT WITHIN A VARIANT CONTEXT FILE===============
    # Tests that the correct donor context is returned
    def test_get_donor_context(self):
        self.assertEqual(self.varcon_file.get_acceptor_context(self.context_id_answer).to_string(),
                         self.acceptor_context_answer.to_string(),
                         "The returned donor context is not what was expected")

    # Tests that None is returned as a non existing donor context is specified
    def test_get_donor_context_none(self):
        self.assertIsNone(self.varcon_file.get_donor_context('22_9411255'),
                          "The requested donor context should not have existed and should have therefore been None")

    # Tests that the correct donor reads are returned
    def test_get_donor_context_reads(self):
        obtained_reads = [x.to_string() for x in self.varcon_file.get_donor_context_reads(self.context_id_answer)]
        answer_reads = [x.to_string for x in self.don_context_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, f"The returned donor context reads should have been "
                             f"{answer_reads}")

    # ===============TESTS FOR READING VARIANT CONTEXT DATA INTO THE VARIANT CONTEXT FILE===============
    # Tests that the variant context file reads a varcon.txt file correctly.
    def test_read_variant_context_file_pos(self):
        read_data_answer = [self.read_varcontxt_entry1, self.read_varcontxt_entry2, self.read_varcontxt_entry3]
        read_data_answer.sort()

        variant_context_file = VariantContextFile(self.read_varcontxt_loc)
        variant_context_file.read_acceptor_context_file(self.acccontxt_loc)
        variant_context_file.read_donor_context_file(self.doncontxt_loc)
        obtained_read_data = [x.to_string() for x in variant_context_file.get_variant_contexts()]
        obtained_read_data.sort()

        print(obtained_read_data)

        #self.assertListEqual(obtained_read_data, read_data_answer,
                             #f"The read variant context data should have been {read_data_answer}")

    # Tests the correct reading of a variant context file with a set sample filter
    def test_read_variant_context_file_samplefilter(self):
        sample_filter = ['SAM001', 'SAM003']
        read_data_answer = [self.read_varcontxt_entry1, self.read_varcontxt_entry3]
        read_data_answer.sort()

        variant_context_file = VariantContextFile(self.read_varcontxt_loc, samplefilter=sample_filter)
        variant_context_file.read_acceptor_context_file(self.acccontxt_loc, samplefilter=sample_filter)
        variant_context_file.read_donor_context_file(self.doncontxt_loc, samplefilter=sample_filter)
        obtained_read_data = [x.to_string() for x in variant_context_file.get_variant_contexts()]
        obtained_read_data.sort()

        self.assertListEqual(obtained_read_data, read_data_answer, "The read variant context data should have been "
                             f"{read_data_answer}")

    # Tests the correct reading of a variant context file with a set context filter
    def test_read_variant_context_file_contextfilter(self):
        context_filter = ['22_9600250', '21_9900000']
        read_data_answer = [self.read_varcontxt_entry2, self.read_varcontxt_entry3]
        read_data_answer.sort()

        variant_context_file = VariantContextFile(self.read_varcontxt_loc, varconfilter=context_filter)
        variant_context_file.read_acceptor_context_file(self.acccontxt_loc, contextfilter=context_filter)
        variant_context_file.read_donor_context_file(self.doncontxt_loc, contextfilter=context_filter)
        obtained_read_data = [x.to_string() for x in variant_context_file.get_variant_contexts()]
        obtained_read_data.sort()

        self.assertListEqual(variant_context_file.get_variant_contexts(), read_data_answer,
                             f"The read variant context data should have been {read_data_answer}")

    # Tests the correct reading of a variant context file with a set chromosome filter
    def test_read_variant_context_file_chromfilter(self):
        chrom_filter = ['22']
        read_data_answer = [self.read_varcontxt_entry2]
        variant_context_file = VariantContextFile(self.read_varcontxt_loc, chromfilter=chrom_filter)
        variant_context_file.read_acceptor_context_file(self.acccontxt_loc, chromfilter=chrom_filter)
        variant_context_file.read_donor_context_file(self.doncontxt_loc, chromfilter=chrom_filter)
        obtained_read_data = [x.to_string() for x in variant_context_file.get_variant_contexts()]
        self.assertListEqual(obtained_read_data, read_data_answer,
                             f"The read variant context data should have been {read_data_answer}")

    #def test_read_variant_context_file_nofile(self):


    #def test_read_variant_context_file_emptyfile(self):

    # ===============TESTS FOR READING ACCEPTOR CONTEXT DATA INTO THE VARIANT CONTEXT FILE===============
    #def test_read_acceptor_context_file(self):
    #def test_read_acceptor_context_file_samplefilter(self):
    #def test_read_acceptor_context_file_contextfilter(self):
    #def test_read_acceptor_context_file_chromfilter(self):
    #def test_read_acceptor_context_file_nofile(self):
    #def test_read_acceptor_context_file_emptyfile(self):

    # ===============TESTS FOR READING DONOR CONTEXT DATA INTO THE VARIANT CONTEXT FILE===============
    #def test_read_donor_context_file(self):
    #def test_read_donor_context_file_samplefilter(self):
    #def test_read_donor_context_file_contextfilter(self):
    #def test_read_donor_context_file_chromfilter(self):
    #def test_read_donor_context_file_nofile(self):
    #def test_read_donor_context_file_emptyfile(self):

    # Tests that a value is indeed in a filter
    def test_passes_filter_pos(self):
        pos_filter_val = 'aap'
        self.assertTrue(self.varcon_file.passes_filter(pos_filter_val, self.filterListToUse),
                        f"The value {pos_filter_val} should have been in the filter list {self.filterListToUse} "
                        "and therefore return True")

    # Tests that a value is not in a filter
    def test_passes_filter_neg(self):
        neg_filter_val = 'jan'
        self.assertFalse(self.varcon_file.passes_filter(neg_filter_val, self.filterListToUse),
                         f"The value {neg_filter_val} should not have been in the filter list {self.filterListToUse} "
                         "and therefore return False")

    # ====================TESTS FOR CHECKING WHETHER A VARIANTS IN A VARIANT CONTEXT====================
    # Tests that variant is in an existing variant context
    def test_variant_is_in_context_pos(self):
        pos_variant_type = 'snp'
        self.assertTrue(self.varcon_file.variant_is_in_context(pos_variant_type, self.context_chrom_answer,
                                                               self.context_origin_answer,
                                                               self.context_origin_answer),
                        f"The variant of type {pos_variant_type} on {self.context_origin_answer}, starting at "
                        f"{self.context_origin_answer} should have been in a context")

    # Tests that a non existing variant is not in a variant context
    def test_variant_is_in_context_neg(self):
        neg_variant_type = 'aap'
        self.assertIsNone(self.varcon_file.variant_is_in_context(neg_variant_type, self.context_chrom_answer,
                                                                 self.context_origin_answer,
                                                                 self.context_origin_answer),
                          f"The variant of type {neg_variant_type} should have returned None")

    def test_snp_variant_is_in_context_pos(self):
        self.assertTrue(self.varcon_file.snp_variant_is_in_context(self.context_chrom_answer,
                                                                   self.context_origin_answer),
                        f"The SNP on chromosome {self.context_chrom_answer} at position {self.context_origin_answer} "
                        "should have been in a variant context")

    def test_snp_variant_is_in_context_neg(self):
        neg_snppos = 325632
        self.assertFalse(self.varcon_file.snp_variant_is_in_context(self.context_chrom_answer, neg_snppos),
                         f"The SNP on chromosome {self.context_chrom_answer} at position {neg_snppos} should not have "
                         "been in any variant context")

    def test_indel_variant_is_in_context_pos(self):
        pos_indel_start = 9411050
        pos_indel_end = 9411150
        self.assertTrue(self.varcon_file.indel_variant_is_in_context(self.context_chrom_answer, pos_indel_start,
                                                                     pos_indel_end),
                        f"The indel on chromosome {self.context_chrom_answer}, starting at {pos_indel_start} and "
                        f"ending at {pos_indel_end} should have been in a variant context")

    def test_indel_variant_is_in_context_neg(self):
        neg_indel_start = 8000000
        neg_indel_end = 8000100
        self.assertFalse(self.varcon_file.indel_variant_is_in_context(self.context_chrom_answer, neg_indel_start,
                                                                      neg_indel_end),
                         f"The indel on chromosome {self.context_chrom_answer}, starting at {neg_indel_start} and "
                         f"ending at {neg_indel_end} should not have been in any variant context")

    def test_context_is_in_variant_context_pos(self):
        variant_context = ["21", 94112080, 9411200, 9411700]
        self.assertTrue(self.varcon_file.context_is_in_variant_context(variant_context),
                        f"Context {variant_context} should been in a variant context")

    def test_context_is_in_variant_context_neg(self):
        variant_context = ["18", 94112080, 9411200, 9411700]
        self.assertFalse(self.varcon_file.context_is_in_variant_context(variant_context),
                         f"Context {variant_context} should not have been in a variant context")

    # ===============TESTS FOF ADDING AND SETTINGS VARIANT, ACCEPTOR AND DONOR CONTEXTS===============
    def test_set_variant_context(self):
        self.varcon_file.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.assertEqual(self.varcon_file.get_variant_context('21_9411000').to_string(),
                         self.setVarContextAnswer.to_string(), "The variant context that was just set and the one "
                                                               "otained are different")

    def test_add_variant_context(self):
        varcon_obj_answer = VariantContext('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100,
                                           self.var_context_a_reads_answer, self.var_context_d_reads_answer,
                                           self.setAccContextAnswer, self.setDonContextAnswer)
        self.varcon_file.add_variant_context('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100,
                                             self.var_context_a_reads_answer, self.var_context_d_reads_answer,
                                             self.setAccContextAnswer, self.setDonContextAnswer)
        self.assertEqual(self.varcon_file.get_variant_context('21_9411000').to_string(),
                         varcon_obj_answer.to_string(), f"The obtained variant context for {self.context_id_answer} "
                         "should have been the same as what was just added")

    def test_set_acceptor_context(self):
        self.varcon_file.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.varcon_file.set_acceptor_context('21_9411000', self.setAccContextAnswer)
        self.assertEqual(self.varcon_file.get_acceptor_context('21_9411000').to_string(),
                         self.setAccContextAnswer.to_string(),
                         "The obtained acceptor context for 21_9411000 should have been the same as what was just set")

    def test_add_acceptor_context(self):
        self.varcon_file.set_variant_context('21_9411000', self.setVarContextAnswer)
        acccon_obj_answer = OverlapContext('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050,
                                           self.acc_context_reads_answer)
        self.varcon_file.add_acceptor_context('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050,
                                              self.acc_context_reads_answer)
        self.assertEqual(self.varcon_file.get_acceptor_context('21_9411000').to_string(),
                         acccon_obj_answer.to_string())

    def test_set_donor_context(self):
        self.varcon_file.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.varcon_file.set_donor_context('21_9411000', self.setDonContextAnswer)
        self.assertEqual(self.varcon_file.get_donor_context('21_9411000').to_string(),
                         self.setDonContextAnswer.to_string(), "The obtained donor context for 21_9411000 should have "
                                                               "been the same as what was just set")

    def test_add_donor_context(self):
        self.varcon_file.set_variant_context('21_9411000', self.setVarContextAnswer)
        doncon_obj_answer = OverlapContext('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100,
                                           self.don_context_reads_answer)
        self.varcon_file.add_donor_context('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100,
                                           self.don_context_reads_answer)
        self.assertEqual(self.varcon_file.get_donor_context('21_9411000').to_string(),
                         doncon_obj_answer.to_string(), "The obtained donor context for '21_9411000' should have been "
                                                        "the same as what was just added")

    # ====================PERFORM THE TESTS FOR SET OPERATIONS ON TWO VARIANT CONTEXT FILES====================
    def test_get_variant_contexts_union(self):
        varcon_union_answer = [self.context_id_answer]
        self.assertListEqual(self.varcon_file.get_variant_contexts_union(self.pos_otherVariantContextFile),
                             varcon_union_answer, f"The variant context union should have been {varcon_union_answer}")

    def get_variant_contexts_intersect_pos(self):
        pos_intersect_answer = [self.context_id_answer]
        self.assertListEqual(self.varcon_file.get_variant_contexts_intersect(self.pos_otherVariantContextFile),
                             pos_intersect_answer, "The variant context intersect should have been "
                             f"{pos_intersect_answer}")

    def test_get_variant_contexts_intersect_neg(self):
        neg_intersect_answer = []
        self.assertListEqual(self.varcon_file.get_variant_contexts_intersect(self.neg_otherVariantContextFile),
                             neg_intersect_answer, f"The variant context intersect should have been empty")

    def test_get_variant_contexts_difference_pos(self):
        pos_difference_answer = ['21_9411259']
        self.assertListEqual(self.varcon_file.get_variant_contexts_difference(self.neg_otherVariantContextFile),
                             pos_difference_answer, "the difference in variant context files should have been "
                             f"{pos_difference_answer} ")

    def test_get_variant_contexts_difference_neg(self):
        neg_difference_answer = []
        self.assertListEqual(self.varcon_file.get_variant_contexts_difference(self.pos_otherVariantContextFile),
                             neg_difference_answer, f"The differences in variant context files should have been empty")

    def test_get_variant_contexts_symmetric_difference_pos(self):
        pos_symdifference_answer = []
        self.assertListEqual(
            self.varcon_file.get_variant_contexts_symmetric_difference(self.pos_otherVariantContextFile),
            pos_symdifference_answer, f"")

    def test_get_variant_contexts_symmetric_difference_neg(self):
        neg_symdifference_answer = [self.context_id_answer, "20_150"]
        self.assertListEqual(self.varcon_file.get_variant_contexts_symmetric_difference(
            self.neg_otherVariantContextFile), neg_symdifference_answer, "The symmetric differences are not the same")

    # ====================PERFORM THE TESTS FOR SETTING UNMAPPED MATE IDS====================
    def test_set_acceptor_context_unmapped_mate_ids(self):
        acu_read_ids = ['acuRead1', 'acuRead3', 'acuRead5']
        self.varcon_file.set_acceptor_context_unmapped_mate_ids(self.context_id_answer, acu_read_ids)
        self.assertListEqual(self.varcon_file.get_acceptor_context_unmapped_mate_ids(self.context_id_answer),
                             acu_read_ids, "The set and returned acceptor context unmapped mate ids should have been "
                             f"{acu_read_ids}")

    def test_set_donor_context_unmapped_mate_ids(self):
        dcu_read_ids = ['dcuRead2', 'dcuRead4', 'dcuRead6']
        self.varcon_file.set_donor_context_unmapped_mate_ids(self.context_id_answer, dcu_read_ids)
        self.assertListEqual(self.varcon_file.get_donor_context_unmapped_mate_ids(self.context_id_answer),
                             dcu_read_ids, "The set and returned donor context unmapped mate ids should have been "
                             f"{dcu_read_ids}")

    def test_set_unmapped_acceptor_mate_ids(self):
        vcau_read_ids = ['vcauRead1', 'vcauRead3', 'vcauRead5']
        self.varcon_file.set_unmapped_acceptor_mate_ids(self.context_id_answer, vcau_read_ids)
        self.assertListEqual(self.varcon_file.get_unmapped_acceptor_mate_ids(self.context_id_answer),
                             vcau_read_ids, "The set and returned unmapped acceptor mate ids should have been: "
                             f"{vcau_read_ids}")

    def test_set_unmapped_donor_mate_ids(self):
        vcdu_read_ids = ['vcduRead2', 'vcduRead4', 'vcduRead6']
        self.varcon_file.set_unmapped_donor_mate_ids(self.context_id_answer, vcdu_read_ids)
        self.assertListEqual(self.varcon_file.get_unmapped_donor_mate_ids(self.context_id_answer),
                             vcdu_read_ids, "The set and returned unmapped donor mate ids should have been: "
                             f"{vcdu_read_ids}")

    # ====================PERFORM THE TESTS FOR WRITING OUTPUT DATA====================
    #def test_write_variant_context_file(self):
    #def test_write_acceptor_context_file(self):
    #def test_write_donor_context_file(self):

    #def test_write_variant_context_stats(self):
    #def test_write_acceptor_context_stats(self):
    #def test_write_donor_context_stats(self):

    #def test_write_left_right_positions(self):
    #def test_write_acceptor_left_right_positions(self):
    #def test_write_donor_left_right_positions(self):

    #def test_write_reads_with_unmapped_mate(self):
    #def test_write_acceptor_unmapped_mate(self):
    #def test_write_donor_unmapped_mate(self):

    # ===============TESTS THE VARIANT CONTEXT COMPARE FUNCTION===============
    #def test_compare(self):
    
    def test_merge_context_windows(self):
        context_window1 = ["21", 150, 100, 200]
        context_window2 = ["21", 200, 150, 250]
        merged_window_answer = ["21", 150, 100, 250]

        self.assertListEqual(self.varcon_file.merge_context_windows(context_window1, context_window2),
                             merged_window_answer,
                             f"The merged context window should have been {merged_window_answer}")

    def test_merge_variant_context_reads(self):
        supply_list = [self.a_read_1, self.a_read_2, self.a_read_2, self.d_read_1, self.d_read_1, self.d_read_2]
        merged_reads_list = [self.a_read_1, self.a_read_2, self.d_read_1, self.d_read_2]

        received_merged_reads = self.varcon_file.merge_variant_context_reads(supply_list)
        merged_reads_answer = "\n"

