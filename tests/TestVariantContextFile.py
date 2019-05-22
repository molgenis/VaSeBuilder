import unittest

# Import the required VaSe classes
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class TestVariantContextFile(unittest.TestCase):
    def setUp(self):
        # Construct the bam reads to use
        self.contextChromAnswer = '21'

        self.vaReadId = 'vaRead1'
        self.vdReadId = 'vdRead1'
        self.aReadId = 'aRead1'
        self.dReadId = 'dRead1'
        self.aReadStartPos_1 = 9411000
        self.aReadStartPos_2 = 9411200
        self.aReadEndPos = 9411350
        self.dReadStartPos_1 = 9411150
        self.dReadStartPos_2 = 9411350
        self.dReadEndPos = 9411500
        self.bamreadLen = 151
        self.aRead1Seq = 'TTTAGATGGG'
        self.aRead2Seq = 'ATTTCTAGTT'
        self.dRead1Seq = 'AGAAAAAGTC'
        self.dRead2Seq = 'TGCCTTTTCA'
        self.aRead1Quals = '=====<===='
        self.aRead2Quals = '>???>?????'
        self.dRead1Quals = '><=???>==<'
        self.dRead2Quals = 'TGCCTTTTCA'
        self.readMapQ = 40
        self.aRead_1 = DonorBamRead(self.aReadId, '1', self.contextChromAnswer, self.aReadStartPos_1, self.bamreadLen, self.aRead1Seq, self.aRead1Quals, self.readMapQ)
        self.aRead_2 = DonorBamRead(self.aReadId, '2', self.contextChromAnswer, self.aReadStartPos_2, self.bamreadLen, self.aRead2Seq, self.aRead2Quals, self.readMapQ)
        self.dRead_1 = DonorBamRead(self.dReadId, '1', self.contextChromAnswer, self.dReadStartPos_1, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        self.dRead_2 = DonorBamRead(self.dReadId, '2', self.contextChromAnswer, self.dReadStartPos_2, self.bamreadLen, self.dRead2Seq, self.dRead2Quals, self.readMapQ)
        self.vaRead_1 = DonorBamRead(self.vaReadId, '1', self.contextChromAnswer, self.aReadStartPos_1, self.bamreadLen, self.aRead1Seq, self.aRead1Quals, self.readMapQ)
        self.vaRead_2 = DonorBamRead(self.vaReadId, '2', self.contextChromAnswer, self.aReadStartPos_2, self.bamreadLen, self.aRead2Seq, self.aRead2Quals, self.readMapQ)
        self.vdRead_1 = DonorBamRead(self.vdReadId, '1', self.contextChromAnswer, self.dReadStartPos_1, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        self.vdRead_2 = DonorBamRead(self.vdReadId, '2', self.contextChromAnswer, self.dReadStartPos_2, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        
        # Create the variables containing context answers and values
        self.contextIdAnswer = '21_9411259'
        self.contextSampleAnswer = 'testsample'
        self.contextOriginAnswer = 9411250
        self.accContextStartAnswer = 9411000
        self.donContextStartAnswer = 9411150
        self.varContextStartAnswer = self.accContextStartAnswer
        self.accContextEndAnswer = 9411350
        self.donContextEndAnswer = 9411500
        self.varContextEndAnswer = self.donContextEndAnswer
        self.accContextReadsAnswer = [self.aRead_1, self.aRead_2]
        self.donContextReadsAnswer = [self.dRead_1, self.dRead_2]
        self.varContextAReadsAnswer = [self.vaRead_1, self.vaRead_2]
        self.varContextDReadsAnswer = [self.vdRead_1, self.vdRead_2]
        self.acceptorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.accContextStartAnswer, self.accContextEndAnswer, self.accContextReadsAnswer)
        self.donorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.donContextStartAnswer, self.donContextEndAnswer, self.donContextReadsAnswer)
        
        # Create the variables containing all the VariantContext answers
        self.variantContextAnswer = VariantContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.varContextStartAnswer, self.varContextEndAnswer, self.varContextAReadsAnswer, self.varContextDReadsAnswer)

        self.setAccContextAnswer = OverlapContext('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050, self.accContextReadsAnswer)
        self.setDonContextAnswer = OverlapContext('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100, self.donContextReadsAnswer)
        self.setVarContextAnswer = VariantContext('21_9411000', 'testanswer', '21', 9411000, 9410900, 9411100, self.varContextAReadsAnswer, self.varContextDReadsAnswer, self.setAccContextAnswer, self.setDonContextAnswer)

        # Create the variables containg the VariantContextFile answers
        self.filterListToUse = ['aap', 'noot', 'mies']
        
        # Create the variables containing the answer of the VariantContextFile
        self.variantContextFile = VariantContextFile()
        self.variantContextFile.set_variant_context(self.contextIdAnswer, self.variantContextAnswer)
        self.variantContextFile.set_acceptor_context(self.contextIdAnswer, self.acceptorContextAnswer)
        self.variantContextFile.set_donor_context(self.contextIdAnswer, self.donorContextAnswer)
        
        # Construct the other VariantContextFile objects to use for the set operation tests
        self.pos_otherVariantContextFile = VariantContextFile()
        self.pos_otherVariantContextFile.set_variant_context(self.contextIdAnswer, self.variantContextAnswer)
        self.neg_otherVariantContextFile = VariantContextFile()
        self.neg_otherVariantContextFile.add_variant_context('20_150', 'testsample2', '20', 150, 0, 350, self.varContextDReadsAnswer, self.varContextAReadsAnswer)
    
    
    
    # ====================PERFORMS THE TESTS FOR THE GETTER METHODS====================
    def test_getVariantContexts(self):
        variant_contexts_answer = [self.variantContextAnswer.to_string()]
        obtained_contexts_answer = [x.to_string() for x in self.variantContextFile.get_variant_contexts()]
        self.assertListEqual(obtained_contexts_answer, variant_contexts_answer, "The returned variant contexts are not what was expected")
    
    def test_getVariantContext(self):
        self.assertEqual(self.variantContextFile.get_variant_context(self.contextIdAnswer).to_string(), self.variantContextAnswer.to_string(), "The returned variant context is not what was expected")
    
    def test_getVariantContext_None(self):
        self.assertIsNone(self.variantContextFile.get_variant_context('22_9411255'), "The requested variant context should not have existed and should have therefore been None")
    
    def test_getAcceptorContext(self):
        self.assertEqual(self.variantContextFile.get_acceptor_context(self.contextIdAnswer).to_string(), self.acceptorContextAnswer.to_string(), "The returned acceptor contexts is not what was expected")
    
    def test_getAcceptorContext_None(self):
        self.assertIsNone(self.variantContextFile.get_acceptor_context('22_9411255'), "The requested acceptor context should not have existed and should have therefore been None")
    
    def test_getDonorContext(self):
        self.assertEqual(self.variantContextFile.get_acceptor_context(self.contextIdAnswer).to_string(), self.acceptorContextAnswer.to_string(), "The returned donor context is not what was expected")
    
    def test_getDonorContext_None(self):
        self.assertIsNone(self.variantContextFile.get_donor_context('22_9411255'), "The requested donor context should not have existed and should have therefore been None")
    
    def test_getAllVariantContextAcceptorReads(self):
        obtained_reads = [x.to_string() for x in self.variantContextFile.get_all_variant_context_acceptor_reads()]
        answer_reads = [x.to_string() for x in self.varContextAReadsAnswer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same data")
    
    def test_getAllVariantContextDonorReads(self):
        obtained_reads = [x.to_string() for x in self.variantContextFile.get_all_variant_context_donor_reads()]
        answer_reads = [x.to_string() for x in self.varContextDReadsAnswer]
        self.assertListEqual(obtained_reads, answer_reads, "Both lists should have contained reads with the exact same data")
    
    def test_getAllVariantContextAcceptorReadIds(self):
        vca_read_ids = [self.vaReadId, self.vaReadId]
        self.assertListEqual(self.variantContextFile.get_all_variant_context_acceptor_read_ids(), vca_read_ids, f"The list of returned acceptor read ids should have been: {vca_read_ids}")

    def test_getAllVariantContextDonorReadIds(self):
        vcd_read_ids = [self.vdReadId, self.vdReadId]
        self.assertListEqual(self.variantContextFile.get_all_variant_context_donor_read_ids(), vcd_read_ids, f"The list of returned donor read ids should have been: {vcd_read_ids}")
    
    
    # ====================PERFORM THE TESTS FOR READING A VARIANT CONTEXT FILE====================
    #def test_readVariantContextFile_pos(self)

    
    def test_passesFilter_pos(self):
        pos_filter_val = 'aap'
        self.assertTrue(self.variantContextFile.passes_filter(pos_filter_val, self.filterListToUse), f"The value {pos_filter_val} should have been in the filter list {self.filterListToUse} and therefore return True")
        
    def test_passesFilter_neg(self):
        neg_filter_val = 'jan'
        self.assertFalse(self.variantContextFile.passes_filter(neg_filter_val, self.filterListToUse), f"The value {neg_filter_val} should not have been in the filter list {self.filterListToUse} and therefore return False")
    
    
    
    # ====================PERFORM THE TESTS FOR IN CONTEXT METHODS====================
    def test_variantIsInContext_pos(self):
        pos_variant_type = 'snp'
        self.assertTrue(self.variantContextFile.variant_is_in_context(pos_variant_type, self.contextChromAnswer, self.contextOriginAnswer, self.contextOriginAnswer), f"The variant of type {pos_variant_type} on {self.contextOriginAnswer}, starting at {self.contextOriginAnswer} should have been in a context")
    
    def test_variantIsInContext_neg(self):
        neg_variant_type = 'aap'
        self.assertIsNone(self.variantContextFile.variant_is_in_context(neg_variant_type, self.contextChromAnswer, self.contextOriginAnswer, self.contextOriginAnswer), f"The variant of type {neg_variant_type} should have returned None")
    
    def test_snpVariantIsInContext_pos(self):
        self.assertTrue(self.variantContextFile.snp_variant_is_in_context(self.contextChromAnswer, self.contextOriginAnswer), f"The SNP on chromosome {self.contextChromAnswer} at position {self.contextOriginAnswer} should have been in a variant context")

    def test_snpVariantIsInContext_neg(self):
        neg_snppos = 325632
        self.assertFalse(self.variantContextFile.snp_variant_is_in_context(self.contextChromAnswer, neg_snppos), f"The SNP on chromosome {self.contextChromAnswer} at position {neg_snppos} should not have been in any variant context")
    
    def test_indelVariantIsInContext_pos(self):
        pos_indel_start = 9411050
        pos_indel_end = 9411150
        self.assertTrue(self.variantContextFile.indel_variant_is_in_context(self.contextChromAnswer, pos_indel_start, pos_indel_end), f"The indel on chromosome {self.contextChromAnswer}, starting at {pos_indel_start} and ending at {pos_indel_end} should have been in a variant context")
    
    def test_indelVariantIsInContext_neg(self):
        neg_indel_start = 8000000
        neg_indel_end = 8000100
        self.assertFalse(self.variantContextFile.indel_variant_is_in_context(self.contextChromAnswer, neg_indel_start, neg_indel_end), f"The indel on chromosome {self.contextChromAnswer}, starting at {neg_indel_start} and ending at {neg_indel_end} should not ahve been in any variant context")
    
    
    
    # ====================PERFORM THE TESTS FOR ADDING CONTEXTS TO THE VARIANT CONTEXT FILE====================
    def test_setVariantContext(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.assertEqual(self.variantContextFile.get_variant_context('21_9411000').to_string(), self.setVarContextAnswer.to_string(), f"The variant context that was just set and the one otained are different")
    
    def test_addVariantContext(self):
        varcon_obj_answer = VariantContext('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100, self.varContextAReadsAnswer, self.varContextDReadsAnswer, self.setAccContextAnswer, self.setDonContextAnswer)
        self.variantContextFile.add_variant_context('21_9411000', 'testanswer', '21', 9411000, 94110900, 9411100, self.varContextAReadsAnswer, self.varContextDReadsAnswer, self.setAccContextAnswer, self.setDonContextAnswer)
        self.assertEqual(self.variantContextFile.get_variant_context('21_9411000').to_string(), varcon_obj_answer.to_string(), f"The obtained variant context for {self.contextIdAnswer} should have been the same as what was just added")
    
    def test_setAcceptorContext(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.variantContextFile.set_acceptor_context('21_9411000', self.setAccContextAnswer)
        self.assertEqual(self.variantContextFile.get_acceptor_context('21_9411000').to_string(), self.setAccContextAnswer.to_string(), f"The obtained acceptor context for 21_9411000 should have been the same as what was just set")
    
    def test_addAcceptorContext(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        acccon_obj_answer = OverlapContext('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050, self.accContextReadsAnswer)
        self.variantContextFile.add_acceptor_context('21_9411000', 'accanswer', '21', 9411000, 94110900, 9411050, self.accContextReadsAnswer)
        self.assertEqual(self.variantContextFile.get_acceptor_context('21_9411000').to_string(), acccon_obj_answer.to_string())
    
    def test_setDonorContext(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        self.variantContextFile.set_donor_context('21_9411000', self.setDonContextAnswer)
        self.assertEqual(self.variantContextFile.get_donor_context('21_9411000').to_string(), self.setDonContextAnswer.to_string(), f"The obtained donor context for 21_9411000 should have been the same as what was just set")
    
    def test_addDonorContext(self):
        self.variantContextFile.set_variant_context('21_9411000', self.setVarContextAnswer)
        doncon_obj_answer = OverlapContext('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100, self.donContextReadsAnswer)
        self.variantContextFile.add_donor_context('21_9411000', 'donanswer', '21', 9411000, 94110950, 9411100, self.donContextReadsAnswer)
        self.assertEqual(self.variantContextFile.get_donor_context('21_9411000').to_string(), doncon_obj_answer.to_string(), f"The obtained donor context for '21_9411000' should have been the same as what was just added")
    
    
    
    # ====================PERFORM THE TESTS FOR SET OPERATIONS ON TWO VARIANT CONTEXT FILES====================
    def test_getVariantContextsUnion(self):
        varcon_union_answer = [self.contextIdAnswer]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_union(self.pos_otherVariantContextFile), varcon_union_answer, f"The variant context union should have been {varcon_union_answer}")
    
    def getVariantContextsIntersect_pos(self):
        pos_intersect_answer = [self.contextIdAnswer]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_intersect(self.pos_otherVariantContextFile), pos_intersect_answer, f"The variant context intersect should have been {pos_intersect_answer}")
    
    def test_getVariantContextsIntersect_neg(self):
        neg_intersect_answer = []
        self.assertListEqual(self.variantContextFile.get_variant_contexts_intersect(self.neg_otherVariantContextFile), neg_intersect_answer, f"The variant context intersect should have been empty")
    
    def test_getVariantContextsDifference_pos(self):
        pos_difference_answer = ['21_9411259']
        self.assertListEqual(self.variantContextFile.get_variant_contexts_difference(self.neg_otherVariantContextFile), pos_difference_answer, f"the difference in variant context files should have been {pos_difference_answer} ")
    
    def test_getVariantContextsDifference_neg(self):
        neg_difference_answer = []
        self.assertListEqual(self.variantContextFile.get_variant_contexts_difference(self.pos_otherVariantContextFile), neg_difference_answer, f"The differences in variant context files should have been empty")
    
    def test_getVariantContextsSymmetricDifference_pos(self):
        pos_symdifference_answer = []
        self.assertListEqual(self.variantContextFile.get_variant_contexts_symmetric_difference(self.pos_otherVariantContextFile), pos_symdifference_answer, f"")
    
    def test_getVariantContextsSymmetricDifference_neg(self):
        neg_symdifference_answer = ['20_150', self.contextIdAnswer]
        self.assertListEqual(self.variantContextFile.get_variant_contexts_symmetric_difference(self.neg_otherVariantContextFile), neg_symdifference_answer, f"The symmetric differences ")



    # ====================PERFORM THE TESTS FOR SETTING UNMAPPED MATE IDS====================
    def test_setAcceptorContextUnmappedMateIds(self):
        acu_read_ids = ['acuRead1', 'acuRead3', 'acuRead5']
        self.variantContextFile.set_acceptor_context_unmapped_mate_ids(self.contextIdAnswer, acu_read_ids)
        self.assertListEqual(self.variantContextFile.get_acceptor_context_unmapped_mate_ids(self.contextIdAnswer), acu_read_ids, f"The set and returned acceptor context unmapped mate ids should have been {acu_read_ids}")

    def test_setDonorContextUnmappedMateIds(self):
        dcu_read_ids = ['dcuRead2', 'dcuRead4', 'dcuRead6']
        self.variantContextFile.set_donor_context_unmapped_mate_ids(self.contextIdAnswer, dcu_read_ids)
        self.assertListEqual(self.variantContextFile.get_donor_context_unmapped_mate_ids(self.contextIdAnswer), dcu_read_ids, f"The set and returned donor context unmapped mate ids should have been {dcu_read_ids}")

    def test_setUnmappedAcceptorMateIds(self):
        vcau_read_ids = ['vcauRead1', 'vcauRead3', 'vcauRead5']
        self.variantContextFile.set_unmapped_acceptor_mate_ids(self.contextIdAnswer, vcau_read_ids)
        self.assertListEqual(self.variantContextFile.get_unmapped_acceptor_mate_ids(self.contextIdAnswer), vcau_read_ids, f"The set and returned unmapped acceptor mate ids should have been: {vcau_read_ids}")

    def test_setUnmappedDonorMateIds(self):
        vcdu_read_ids = ['vcduRead2', 'vcduRead4', 'vcduRead6']
        self.variantContextFile.set_unmapped_donor_mate_ids(self.contextIdAnswer, vcdu_read_ids)
        self.assertListEqual(self.variantContextFile.get_unmapped_donor_mate_ids(self.contextIdAnswer), vcdu_read_ids, f"The set and returned unmapped donor mate ids should have been: {vcdu_read_ids}")



    # ====================PERFORM THE TESTS FOR WRITING OUTPUT DATA====================
