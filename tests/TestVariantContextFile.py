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

        # Create the variables containg the VariantContextFile answers
        self.filterListToUse = ["aap", "noot", "mies"]
        
        # Create the variables containing the answer of the VariantContextFile
        self.variantContextFile = VariantContextFile()
        self.variantContextFile.set_variant_context(self.context_id_answer, self.variantContextAnswer)
        self.variantContextFile.set_acceptor_context(self.context_id_answer, self.acceptor_context_answer)
        self.variantContextFile.set_donor_context(self.context_id_answer, self.donor_context_answer)
        
        # Construct the other VariantContextFile objects to use for the set operation tests
        self.pos_otherVariantContextFile = VariantContextFile()
        self.pos_otherVariantContextFile.set_variant_context(self.context_id_answer, self.variantContextAnswer)
        self.neg_otherVariantContextFile = VariantContextFile()
        self.neg_otherVariantContextFile.add_variant_context('20_150', 'testsample2', '20', 150, 0, 350,
                                                             self.var_context_d_reads_answer,
                                                             self.var_context_a_reads_answer)

    # ====================PERFORMS THE TESTS FOR THE GETTER METHODS====================
    def test_get_variant_contexts(self):
        variant_contexts_answer = [self.variantContextAnswer.to_string()]
        obtained_contexts_answer = [x.to_string() for x in self.variantContextFile.get_variant_contexts()]
        self.assertListEqual(obtained_contexts_answer, variant_contexts_answer,
                             "The returned variant contexts are not what was expected")

    def test_get_variant_context(self):
        self.assertEqual(self.variantContextFile.get_variant_context(self.context_id_answer).to_string(),
                         self.variantContextAnswer.to_string(), "The returned variant context is not what was expected")

    def test_get_variant_context_none(self):
        self.assertIsNone(self.variantContextFile.get_variant_context('22_9411255'),
                          "The requested variant context should not have existed and should have therefore been None")

    def test_get_acceptor_context(self):
        self.assertEqual(self.variantContextFile.get_acceptor_context(self.context_id_answer).to_string(),
                         self.acceptor_context_answer.to_string(),
                         "The returned acceptor contexts is not what was expected")

    def test_get_acceptor_context_none(self):
        self.assertIsNone(self.variantContextFile.get_acceptor_context("22_9411255"),
                          "The requested acceptor context should not have existed and should have therefore been None")

    def test_get_donor_context(self):
        self.assertEqual(self.variantContextFile.get_acceptor_context(self.context_id_answer).to_string(),
                         self.acceptor_context_answer.to_string(),
                         "The returned donor context is not what was expected")

    def test_get_donor_context_none(self):
        self.assertIsNone(self.variantContextFile.get_donor_context('22_9411255'),
                          "The requested donor context should not have existed and should have therefore been None")

    def test_get_all_variant_context_acceptor_reads(self):
        obtained_reads = [x.to_string() for x in self.variantContextFile.get_all_variant_context_acceptor_reads()]
        answer_reads = [x.to_string() for x in self.var_context_a_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same "
                                                           "data")

    def test_get_all_variant_context_donor_reads(self):
        obtained_reads = [x.to_string() for x in self.variantContextFile.get_all_variant_context_donor_reads()]
        answer_reads = [x.to_string() for x in self.var_context_d_reads_answer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same "
                                                           "data")

    def test_get_all_variant_context_acceptor_read_ids(self):
        vca_read_ids = [self.va_read_id]
        self.assertListEqual(self.variantContextFile.get_all_variant_context_acceptor_read_ids(), vca_read_ids,
                             f"The list of returned acceptor read ids should have been: {vca_read_ids}")

    def test_get_all_variant_context_donor_read_ids(self):
        vcd_read_ids = [self.vd_read_id]
        self.assertListEqual(self.variantContextFile.get_all_variant_context_donor_read_ids(), vcd_read_ids,
                             f"The list of returned donor read ids should have been: {vcd_read_ids}")

    # ====================PERFORM THE TESTS FOR READING A VARIANT CONTEXT FILE====================
    #def test_readVariantContextFile_pos(self)

    def test_passes_filter_pos(self):
        pos_filter_val = 'aap'
        self.assertTrue(self.variantContextFile.passes_filter(pos_filter_val, self.filterListToUse),
                        f"The value {pos_filter_val} should have been in the filter list {self.filterListToUse} "
                        "and therefore return True")

    def test_passes_filter_neg(self):
        neg_filter_val = 'jan'
        self.assertFalse(self.variantContextFile.passes_filter(neg_filter_val, self.filterListToUse),
                         f"The value {neg_filter_val} should not have been in the filter list {self.filterListToUse} "
                         "and therefore return False")

    # ====================PERFORM THE TESTS FOR IN CONTEXT METHODS====================
    def test_variant_is_in_context_pos(self):
        pos_variant_type = 'snp'
        self.assertTrue(self.variantContextFile.variant_is_in_context(pos_variant_type, self.context_chrom_answer,
                                                                      self.context_origin_answer,
                                                                      self.context_origin_answer),
                        f"The variant of type {pos_variant_type} on {self.context_origin_answer}, starting at "
                        f"{self.context_origin_answer} should have been in a context")

    def test_variant_is_in_context_neg(self):
        neg_variant_type = 'aap'
        self.assertIsNone(self.variantContextFile.variant_is_in_context(neg_variant_type, self.context_chrom_answer,
                                                                        self.context_origin_answer,
                                                                        self.context_origin_answer),
                          f"The variant of type {neg_variant_type} should have returned None")

    def test_snp_variant_is_in_context_pos(self):
        self.assertTrue(self.variantContextFile.snp_variant_is_in_context(self.context_chrom_answer,
                                                                          self.context_origin_answer),
                        f"The SNP on chromosome {self.context_chrom_answer} at position {self.context_origin_answer} "
                        "should have been in a variant context")

    def test_snp_variant_is_in_context_neg(self):
        neg_snppos = 325632
        self.assertFalse(self.variantContextFile.snp_variant_is_in_context(self.context_chrom_answer, neg_snppos),
                         f"The SNP on chromosome {self.context_chrom_answer} at position {neg_snppos} should not have "
                         "been in any variant context")

    def test_indel_variant_is_in_context_pos(self):
        pos_indel_start = 9411050
        pos_indel_end = 9411150
        self.assertTrue(self.variantContextFile.indel_variant_is_in_context(self.context_chrom_answer, pos_indel_start,
                                                                            pos_indel_end),
                        f"The indel on chromosome {self.context_chrom_answer}, starting at {pos_indel_start} and "
                        f"ending at {pos_indel_end} should have been in a variant context")

    def test_indel_variant_is_in_context_neg(self):
        neg_indel_start = 8000000
        neg_indel_end = 8000100
        self.assertFalse(self.variantContextFile.indel_variant_is_in_context(self.context_chrom_answer, neg_indel_start,
                                                                             neg_indel_end),
                         f"The indel on chromosome {self.context_chrom_answer}, starting at {neg_indel_start} and "
                         f"ending at {neg_indel_end} should not ahve been in any variant context")

    # ====================PERFORM THE TESTS FOR ADDING CONTEXTS TO THE VARIANT CONTEXT FILE====================
    def test_set_variant_context(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.assertEqual(self.variantContextFile.get_variant_context('21_9411000').to_string(),
                         self.setVarContextAnswer.to_string(), "The variant context that was just set and the one "
                                                               "otained are different")

    def test_add_variant_context(self):
        varcon_obj_answer = VariantContext('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100,
                                           self.var_context_a_reads_answer, self.var_context_d_reads_answer,
                                           self.setAccContextAnswer, self.setDonContextAnswer)
        self.variantContextFile.add_variant_context('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100,
                                                    self.var_context_a_reads_answer, self.var_context_d_reads_answer,
                                                    self.setAccContextAnswer, self.setDonContextAnswer)
        self.assertEqual(self.variantContextFile.get_variant_context('21_9411000').to_string(),
                         varcon_obj_answer.to_string(), f"The obtained variant context for {self.context_id_answer} "
                         "should have been the same as what was just added")

    def test_set_acceptor_context(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.variantContextFile.set_acceptor_context('21_9411000', self.setAccContextAnswer)
        self.assertEqual(self.variantContextFile.get_acceptor_context('21_9411000').to_string(),
                         self.setAccContextAnswer.to_string(),
                         "The obtained acceptor context for 21_9411000 should have been the same as what was just set")

    def test_add_acceptor_context(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        acccon_obj_answer = OverlapContext('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050,
                                           self.acc_context_reads_answer)
        self.variantContextFile.add_acceptor_context('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050,
                                                     self.acc_context_reads_answer)
        self.assertEqual(self.variantContextFile.get_acceptor_context('21_9411000').to_string(),
                         acccon_obj_answer.to_string())

    def test_set_donor_context(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.variantContextFile.set_donor_context('21_9411000', self.setDonContextAnswer)
        self.assertEqual(self.variantContextFile.get_donor_context('21_9411000').to_string(),
                         self.setDonContextAnswer.to_string(), "The obtained donor context for 21_9411000 should have "
                                                               "been the same as what was just set")

    def test_add_donor_context(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        doncon_obj_answer = OverlapContext('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100,
                                           self.don_context_reads_answer)
        self.variantContextFile.add_donor_context('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100,
                                                  self.don_context_reads_answer)
        self.assertEqual(self.variantContextFile.get_donor_context('21_9411000').to_string(),
                         doncon_obj_answer.to_string(), "The obtained donor context for '21_9411000' should have been "
                                                        "the same as what was just added")

    # ====================PERFORM THE TESTS FOR SET OPERATIONS ON TWO VARIANT CONTEXT FILES====================
    def test_get_variant_contexts_union(self):
        varcon_union_answer = [self.context_id_answer]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_union(self.pos_otherVariantContextFile),
                             varcon_union_answer, f"The variant context union should have been {varcon_union_answer}")

    def get_variant_contexts_intersect_pos(self):
        pos_intersect_answer = [self.context_id_answer]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_intersect(self.pos_otherVariantContextFile),
                             pos_intersect_answer, "The variant context intersect should have been "
                             f"{pos_intersect_answer}")

    def test_get_variant_contexts_intersect_neg(self):
        neg_intersect_answer = []
        self.assertListEqual(self.variantContextFile.get_variant_contexts_intersect(self.neg_otherVariantContextFile),
                             neg_intersect_answer, f"The variant context intersect should have been empty")

    def test_get_variant_contexts_difference_pos(self):
        pos_difference_answer = ['21_9411259']
        self.assertListEqual(self.variantContextFile.get_variant_contexts_difference(self.neg_otherVariantContextFile),
                             pos_difference_answer, "the difference in variant context files should have been "
                             f"{pos_difference_answer} ")

    def test_get_variant_contexts_difference_neg(self):
        neg_difference_answer = []
        self.assertListEqual(self.variantContextFile.get_variant_contexts_difference(self.pos_otherVariantContextFile),
                             neg_difference_answer, f"The differences in variant context files should have been empty")

    def test_get_variant_contexts_symmetric_difference_pos(self):
        pos_symdifference_answer = []
        self.assertListEqual(
            self.variantContextFile.get_variant_contexts_symmetric_difference(self.pos_otherVariantContextFile),
            pos_symdifference_answer, f"")

    def test_get_variant_contexts_symmetric_difference_neg(self):
        neg_symdifference_answer = [self.context_id_answer, "20_150"]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_symmetric_difference(
            self.neg_otherVariantContextFile), neg_symdifference_answer, "The symmetric differences are not the same")

    # ====================PERFORM THE TESTS FOR SETTING UNMAPPED MATE IDS====================
    def test_set_acceptor_context_unmapped_mate_ids(self):
        acu_read_ids = ['acuRead1', 'acuRead3', 'acuRead5']
        self.variantContextFile.set_acceptor_context_unmapped_mate_ids(self.context_id_answer, acu_read_ids)
        self.assertListEqual(self.variantContextFile.get_acceptor_context_unmapped_mate_ids(self.context_id_answer),
                             acu_read_ids, "The set and returned acceptor context unmapped mate ids should have been "
                             f"{acu_read_ids}")

    def test_set_donor_context_unmapped_mate_ids(self):
        dcu_read_ids = ['dcuRead2', 'dcuRead4', 'dcuRead6']
        self.variantContextFile.set_donor_context_unmapped_mate_ids(self.context_id_answer, dcu_read_ids)
        self.assertListEqual(self.variantContextFile.get_donor_context_unmapped_mate_ids(self.context_id_answer),
                             dcu_read_ids, "The set and returned donor context unmapped mate ids should have been "
                             f"{dcu_read_ids}")

    def test_set_unmapped_acceptor_mate_ids(self):
        vcau_read_ids = ['vcauRead1', 'vcauRead3', 'vcauRead5']
        self.variantContextFile.set_unmapped_acceptor_mate_ids(self.context_id_answer, vcau_read_ids)
        self.assertListEqual(self.variantContextFile.get_unmapped_acceptor_mate_ids(self.context_id_answer),
                             vcau_read_ids, "The set and returned unmapped acceptor mate ids should have been: "
                             f"{vcau_read_ids}")

    def test_set_unmapped_donor_mate_ids(self):
        vcdu_read_ids = ['vcduRead2', 'vcduRead4', 'vcduRead6']
        self.variantContextFile.set_unmapped_donor_mate_ids(self.context_id_answer, vcdu_read_ids)
        self.assertListEqual(self.variantContextFile.get_unmapped_donor_mate_ids(self.context_id_answer),
                             vcdu_read_ids, "The set and returned unmapped donor mate ids should have been: "
                             f"{vcdu_read_ids}")

    # ====================PERFORM THE TESTS FOR WRITING OUTPUT DATA====================
