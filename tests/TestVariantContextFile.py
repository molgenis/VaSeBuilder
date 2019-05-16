import unittest

# Import the required VaSe classes
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class TestVariantContextFile(unittest.TestCase):
    def setUp(self):
        # Construct the bam reads to use
        self.vaReadId = 'vaRead1'
        self.vdReadId = 'vdRead1'
        self.aReadId = 'aRead1'
        self.dReadId = 'dRead1'
        self.aReadStartPos_1 = 9411000
        self.aReadStartPos_2 = 9411200
        self.aReadEndPos = 9411350
        self.dReadStartPos_1 = 9411150
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
        self.aRead_1 = DonorBamRead(self.aReadId '1', self.contextChromAnswer, self.aReadStartPos_1, self.bamreadLen, self.aRead1Seq, self.aRead1Quals, self.readMapQ)
        self.aRead_2 = DonorBamRead(self.aReadId, '2', self.contextChromAnswer, self.aReadStartPos_2, self.bamreadLen, self.aRead2Seq, self.aRead2Quals, self.readMapQ)
        self.dRead_1 = DonorBamRead(self.dReadId '1', self.contextChromAnswer, self.dReadStartPos_1, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        self.dRead_2 = DonorBamRead(self.dReadId, '2', self.contextChromAnswer, self.dReadStartPos_2, self.bamreadLen, self.dRead2Seq, self.dRead2Quals, self.readMapQ)
        self.vaRead_1 = DonorBamRead(self.vaReadId, '1', self.contextChromAnswer, self.aReadStartPos_1, self.bamreadLen, self.aRead1Seq, self.aRead1Quals, self.readMapQ)
        self.vaRead_2 = DonorBamRead(self.vaReadId, '2', self.contextChromAnswer, self.aReadStartPos_2, self.bamreadLen, self.aRead2Seq, self.aRead2Quals, self.readMapQ)
        self.vdRead_1 = DonorBamRead(self.vdReadId '1', self.contextChromAnswer, self.dReadStartPos_1, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        self.vdRead_2 = DonorBamRead(self.vdReadId '2', self.contextChromAnswer, self.dReadStartPos_2, self.bamreadLen, self.dRead1Seq, self.dRead1Quals, self.readMapQ)
        
        # Create the variables containing context answers and values
        self.contextIdAnswer = '21_9411259'
        self.contextSampleAnswer = 'testsample'
        self.contextChromAnswer = 21
        self.contextOriginAnswer = 9411250
        self.accContextStartAnswer = 9411000
        self.donContextStartAnswer = 9411150
        self.varContextStartAnswer = self.accContextStartAnswer
        self.accContextEndAnswer = 9411350
        self.donContextEndAnswer = 9411500
        self.varContextEndAnswer = self.donContextEndAnswer
        self.accContextReadsAnswer = [self.aRead1, self.aRead2]
        self.donContextReadsAnswer = [self.dRead1, self.dRead2]
        self.varContextAReadsAnswer = [self.vaRead_1, self.vaRead_2]
        self.varContextDReadsAnswer = [self.vdRead_1, self.vdRead_2]
        self.acceptorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.accContextStartAnswer, self.accContextEndAnswer, self.accContextReadsAnswer)
        self.donorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.donContextStartAnswer, self.donContextEndAnswer, self.donContextReadsAnswer)
        
        # Create the variables containing all the VariantContext answers
        self.variantContextAnswer = VariantContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.varContextStartAnswer, self.varContextEndAnswer, self.varContextAReadsAnswer, self.varContextDReadsAnswer)
        self.pos_varconAnswer = VariantContext()
        self.neg_varconAnswer = VariantContext()
        
        # Create the variables containg the VariantContextFile answers
        self.variantContextsAnswer = {self.contextIdAnswer:self.variantContextAnswer}
        self.posFilterValToUse = 'aap'
        self.negFilterValToUse = 'jan'
        self.filterListToUse = ['aap', 'noot', 'mies']
        
        # Create the variables containing the answer of the VariantContextFile
        self.variantContextFile = VariantContextFile()
        
        # Construct the other VariantContextFile objects to use for the set operation tests
        self.pos_otherVariantContextFile = VariantContextFile()
        self.pos_otherVariantContextFile.setVariantContext(self.contextIdAnswer, self.variantContextAnswer)
        self.neg_otherVariantContextFile = VariantContextFile()
        self.neg_otherVariantContextFile.addVariantContext('20_150', 'testsample2', '20', 150, 0, 350, self.varContextDReadsAnswer, self.varContextAReadsAnswer)
    
    
    
    # ====================PERFORMS THE TESTS FOR THE GETTER METHODS====================
    def test_getVariantContexts(self):
        self.assertDictEqual(self.variantContextFile.getVariantContexts(), self.variantContextsAnswer, "The returned variant contexts are not what was expected")
    
    def test_getVariantContext(self):
        self.assertEqual(self.variantContextFile.getVariantContext(self.contextIdAnswer).toString(), self.variantContextAnswer.toString(), "The returned variant context is not what was expected")
    
    def test_getVariantContext_None(self):
        self.assertIsNone(self.variantContextFile.getVariantContext('22_9411255'), "The requested variant context should not have existed and should have therefore been None")
    
    def test_getAcceptorContext(self):
        self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer).toString(), self.acceptorContextAnswer.toString(), "The returned acceptor contexts is not what was expected")
    
    def test_getAcceptorContext_None(self):
        self.assertIsNone(self.variantContextFile.getAcceptorContext('22_9411255'), "The requested acceptor context should not have existed and should have therefore been None")
    
    def test_getDonorContext(self):
        self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer), self.donorContextAnswer.toString(), "The returned donor context is not what was expected")
    
    def test_getDonorContext_None(self):
        self.assertIsNone(self.variantContextFile.getDonorContext('22_9411255'), "The requested donor context should not have existed and should have therefore been None")
    
    def test_getAllVariantContextAcceptorReads(self):
        self.assertListEqual(self.variantContextFile.getAllVariantContextAcceptorReads(), self.allVarconAcceptorReadsAnswer, "")
    
    def test_getAllVariantContextDonorReads(self):
        self.assertListEqual(self.variantContextFile.getAllVariantContextAcceptorReads(), self.allVarconDonorReadsAnswer, "")
    
    def test_getAllVariantContextAcceptorReadIds(self):
    
    
    def test_getAllVariantContextDonorReadIds(self):
        
    
    
    # ====================PERFORM THE TESTS FOR READING A VARIANT CONTEXT FILE====================
    #def test_readVariantContextFile_pos(self)
    
    def test_readVariantContextFile_neg(self):
        with self.assertRaises(IOError) as cm:
            self.variantContextFile.readVariantContextFile('nonexisting.txt')
        raisedException = cm.exception
        self.assertEqual(raisedException.error_code, 3)
    
    def test_passesFilter_pos(self):
        posFilterValToUse = 'aap'
        self.assertTrue(self.variantContextFile.passesFilter(posFilterValToUse, self.filterListToUse), f"The value {posFilterValToUse} should have been in the filter list {self.filterListToUse} and therefore return True")
        
    def test_passesFilter_neg(self):
        negFilterValToUse = 'jan'
        self.assertFalse(self.variantContextFile.passesFilter(negFilterValToUse, self.filterListToUse), f"The value {negFilterValToUse} should not have been in the filter list {self.filterListToUse} and therefore return False")
    
    
    
    # ====================PERFORM THE TESTS FOR IN CONTEXT METHODS====================
    def test_variantIsInContext_pos(self):
        pos_variantTypeToUse = 'snp'
        self.assertTrue(self.variantContextFile.variantIsInContext(pos_variantTypeToUse, self.variantChromToUse, self.pos_snpPosToUse, self.pos_snpPosToUse), f"The variant of type {pos_variantTypeToUse} on {self.variantChromToUse}, starting at {self.pos_snpPosToUse} should have been in a context")
    
    def test_variantIsInContext_neg(self):
        neg_variantTypeToUse = 'aap'
        self.assertIsNone(self.variantContextFile.variantIsInContext(neg_variantTypeToUse, self.variantChromToUse, self.pos_snpPosToUse, self.pos_snpPosToUse), f"The variant of type {neg_variantTypeToUse} should have returned None")
    
    def test_snpVariantIsInContext_pos(self):
        self.assertTrue(self.variantContextFile.snpVariantIsInContext(self.variantChromToUse, self.pos_snpPosToUse), f"The SNP on chromosome {self.variantChromToUse} at position {self.pos_snpPosToUse} should have been in a variant context")
    
    def test_snpVariantIsInContext_neg(self):
        self.assertFalse(self.variantContextFile.snpVariantIsInContext(self.variantChromToUse, self.neg_snpPosToUse), f"The SNP on chromosome {self.variantChromToUse} at position {self.neg_snpPosToUse} should not have been in any variant context")
    
    def test_indelVariantIsInContext_pos(self):
        pos_indelStartPosToUse = 0
        pos_indelEndPosToUse = 0
        self.assertTrue(self.variantContextFile.indelVariantIsInContext(self.variantChromToUse, pos_indelStartPosToUse, pos_indelEndPosToUse), f"The indel on chromosome {self.variantChromToUse}, starting at {pos_indelStartPosToUse} and ending at {pos_indelEndPosToUse} should have been in a variant context")
    
    def test_indelVariantIsInContext_neg(self):
        neg_indelStartPosToUse = 8000000
        neg_indelEndPosToUse = 8000100
        self.assertFalse(self.variantContextFile.indelVariantIsInContext(self.variantChromToUser, neg_indelStartPosToUse, neg_indelEndPosToUse), f"The indel on chromosome {self.variantChromToUse}, starting at {neg_indelStartPosToUse} and ending at {neg_indelEndPosToUse} should not ahve been in any variant context")
    
    
    
    # ====================PERFORM THE TESTS FOR ADDING CONTEXTS TO THE VARIANT CONTEXT FILE====================
    def test_setVariantContext(self):
        self.variantContextFile.setVariantContext(self.addVarContextObjAnswer)
        self.assertEqual(self.variantContextFile.getVariantContext().toString(), varcon.toString(), "")
    
    def test_addVariantContext(self):
        self.variantContextFile.addVariantContext(self.addContextIdAnswer, self.addContextChromAnswer, self.addContextOriginAnswer, self.addVarContextStartAnswer, self.addVarContextEndAnswer, self.addContextReadsAnswer, self.addContextReadsAnswer)
        self.assertEqual(self.addVarContextObjAnswer.toString(), self.addVarContextObjAnswer.toString(), f"The obtained variant context for {self.addContextIdAnswer} should have been the same as what was just added")
    
    def test_setAcceptorContext(self):
        self.variantContextFile(self.contextIdAnswer, self.addAccContextObjAnswer)
        self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer).toString(), self.addAccContextObjAnswer.toString(), f"The obtained acceptor context for {self.contextIdAnswer} should have been the same as what was just set")
    
    def test_addAcceptorContext(self):
        self.variantContextFile.addAcceptorContext(self.addContextIdAnswer, self.addContextChromAnswer, self.addContextOriginAnswer, self.addAccContextStartAnswer, self.addAccContextEndAnswer)
        self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer), self.addAccContextObjAnswer.toString())
    
    def test_setDonorContext(self):
        self.variantContextFile.setDonorContext(self.contextIdAnswer, donCon)
        self.assertEqual(self.variantContextFile.getDonorContext(self.contextIdAnswer).toString(), self.addDonContextObjAnswer.toString(), f"The obtained donor context for {self.contextIdAnswer} should have been the same as what was just set")
    
    def test_addDonorContext(self):
        self.variantContextFile.addDonorContext(self.addContextIdAnswer, self.addContextChromAnswer, self.addContextOriginAnswer, self.addDonContextStartAnswer, self.addDonContextEndAnswer)
        self.assertEqual(self.variantContextFile.getDonorContext(self.contextIdAnswer).toString(), self.addDonContextObjAnswer.toString(), f"The obtained donor context for {self.contextIdAnswer} should have been the same as what was just added")
    
    
    
    # ====================PERFORM THE TESTS FOR SET OPERATIONS ON TWO VARIANT CONTEXT FILES====================
    def test_getVariantContextsUnion(self, otherVarconFile):
        varconUnionAnswer = [self.contextIdAnswer]
        self.assertListEqual(self.variantContextFile.getVariantContextsUnion(, self.pos_otherVariantContextFile), varconUnionAnswer, f"The variant context union should have been {pos_varconUnionAnswer}")
    
    def getVariantContextsIntersect_pos(self, otherVarconFile):
        pos_varconIntersectAnswer = [self.contextIdAnswer]
        self.asserListEqual(self.variantContextFile().getVariantContextsIntersect(self.pos_otherVariantContextFile), pos_varconIntersectAnswer, f"The variant context intersect should have been {pos_varconIntersectAnswer}")
    
    def test_getVariantContextsIntersect_neg(self, otherVarconFile):
        neg_varconIntersectAnswer = []
        self.assertListEqual(self.variantContextFile.getVariantContextsIntersect(self._neg_otherVariantContextFile), neg_varconIntersectAnswer, f"The variant context intersect should have been empty")
    
    def getVariantContextsDifference_pos(self, otherVarconFile):
        pos_varconDiffAnswer = ['20_150']
        self.assertListEqual(self.variantContextFile.getVariantContextsDifference(), pos_varconDiffAnswer, f"the difference in variant context files should have been {pos_varconDiffAnswer} ")
    
    def test_getVariantContextsDifference_neg(self, otherVarconFile):
        neg_varconDiffAnswer = []
        self.asserListEqual(self.variantContextFile.getVariantContextsDifference(), neg_varconDiffAnswer, f"The differences in variant context files should have been empty")
    
    def getVariantContextsSymmetricDifference_pos(self, otherVarconFile):
        pos_varconSymDiffAnswer = []
        self.assertListEqual(self.variantContextFile.getVariantContextsSymmetricDifference(), pos_varconSymDiffAnswer, f"")
    
    def getVariantContextsSymmetricDifference_neg(self, otherVarconFile):
        neg_varconSymDiffAnswer = [self.contextIdAnswer, '20_150']
        self.assertListEqual(self.variantContextFile.getVariantContextsSymmetricDifference(), neg_varconSymDiffAnswer, f"The symmetric differences ")
