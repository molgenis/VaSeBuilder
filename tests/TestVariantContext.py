import unittest
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext

class TestVariantContext(unittest.TestCase):
    def setUp(self):
        # Create the variables saving the context BAM reads
        self.readIdAnswer = 'HHKY2CCXX160108:1:2122:24160:2522'
        self.readDonorIdAnswer = 'HHKY2CCXX160108:1:2122:24160:2555'
        self.readPnAnswer = '1'
        self.readDonorPnAnswer = '2'
        self.readChromAnswer = '21'
        self.readPosAnswer = 9411193
        self.readLenAnswer = 151
        self.readDonorLenAnswer = 108
        self.readSeqAnswer = 'AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGACTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA'
        self.readQualsAnswer = '><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>>==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>'
        self.readMapQAnswer = 40
        self.acceptorRead = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
        self.donorRead = DonorBamRead(self.readDonorIdAnswer, self.readDonorPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readDonorLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
        
        self.acceptorReadIdsAnswer = ['HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522']
        self.donorReadIdsAnswer = ['HHKY2CCXX160108:1:2122:24160:2555', 'HHKY2CCXX160108:1:2122:24160:2555', 'HHKY2CCXX160108:1:2122:24160:2555']
        self.readLensAnswer = [self.readLenAnswer, self.readLenAnswer, self.readLenAnswer]
        self.donorReadLensAnswer = [self.readDonorLenAnswer, self.readDonorLenAnswer, self.readDonorLenAnswer]
        
        # Create the variables containing the context answers
        self.contextIdAnswer = '21_9411259'
        self.contextSampleAnswer = 'testsample'
        self.contextChromAnswer = '21'
        self.contextOriginAnswer = 9411259
        self.acceptorStartAnswer = 9411210
        self.donorStartAnswer =  9411193
        self.contextStartAnswer = 9411193
        self.acceptorEndAnswer = 9411344
        self.donorEndAnswer = 9411301
        self.contextEndAnswer = 9411344
        self.avgMedLenAnswer = []
        self.avgMedQualAnswer = []
        self.avgMedMapQAnswer = []
        
        # Creat the variables containing the acceptor context answers
        self.acceptorReadsAnswer = [self.acceptorRead, self.acceptorRead, self.acceptorRead]
        self.acceptorReadStartsAnswer = [self.readPosAnswer, self.readPosAnswer, self.readPosAnswer]
        self.acceptorReadEndsAnswer = [9411344, 9411344, 9411344]
        self.acceptorRightPosAnswer = []
        
        # Create the variables containing the donor context answers
        self.donorReadsAnswer = [self.donorRead, self.donorRead, self.donorRead]
        self.donorReadStartsAnswer = [9411193, 9411193, 9411193]
        self.donorReadEndsAnswer = [9411301, 9411301, 9411301]
        self.donorLeftPosAnswer = []
        
        # Create the variables containing the accept, donor and variant context
        self.acceptorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.acceptorStartAnswer, self.acceptorEndAnswer, self.acceptorReadsAnswer)
        self.donorContextAnswer = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.donorStartAnswer, self.donorEndAnswer, self.donorReadsAnswer)
        self.variantContext = VariantContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.contextStartAnswer, self.contextEndAnswer, self.acceptorReadsAnswer, self.donorReadsAnswer, self.acceptorContextAnswer, self.donorContextAnswer)
        
        # Create the variables containing the string output answers
        self.toStatisticsAnswer = f"{self.contextIdAnswer}\t"
    
    
    
    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_getVariantContextId(self):
        self.assertEqual(self.variantContext.get_variant_context_id(), self.contextIdAnswer, f"The variant context id should have been {self.contextIdAnswer}")
    
    def test_getVariantContextSample(self):
        self.assertEqual(self.variantContext.get_variant_context_sample(), self.contextSampleAnswer, f"The variant context sample should have been {self.contextSampleAnswer}")
    
    def test_getVariantContextChrom(self):
        self.assertEqual(self.variantContext.get_variant_context_chrom(), self.contextChromAnswer, f"The variant context chromosome should have been {self.contextChromAnswer}")
    
    def test_getVariantContextOrigin(self):
        self.assertEqual(self.variantContext.get_variant_context_origin(), self.contextOriginAnswer, f"The variant context origin position should have been {self.contextOriginAnswer}")
    
    def test_getVariantContextStart(self):
        self.assertEqual(self.variantContext.get_variant_context_start(), self.contextStartAnswer, f"The variant context start position should have been {self.contextStartAnswer}")
    
    def test_getVariantContextEnd(self):
        self.assertEqual(self.variantContext.get_variant_context_end(), self.contextEndAnswer, f"The variant context end position should have been {self.contextEndAnswer}")
    
    def test_getVariantContextAcceptorReads(self):
        self.assertListEqual(self.variantContext.get_acceptor_reads(), self.acceptorReadsAnswer, f"The variant context acceptor reads are not what was expected")
    
    def test_getVariantContextDonorReads(self):
        self.assertListEqual(self.variantContext.get_donor_reads(), self.donorReadsAnswer, f"The variant context donor reads are not what was expected")
    
    def test_getAcceptorContext(self):
        self.assertEqual(self.variantContext.get_acceptor_context(), self.acceptorContextAnswer, f"The acceptor context is not what was expected")
    
    def test_getDonorContext(self):
        self.assertEqual(self.variantContext.get_donor_context(), self.donorContextAnswer, f"The donor context is not what was expected")
    
    def test_getUnmappedAcceptorMateIds(self):
        vuaReadIds = ['vuaRead1', 'vuaRead2', 'vuaRead3']
        self.variantContext.set_unmapped_acceptor_mate_ids(vuaReadIds)
        self.assertListEqual(self.variantContext.get_unmapped_acceptor_mate_ids(), vuaReadIds, f"The variant context unmapped acceptor mate ids should have been {vuaReadIds}")
    
    def test_getUnmappedDonorMateIds(self):
        vudReadIds = ['vudRead1', 'vudRead2', 'vudRead3']
        self.variantContext.set_unmapped_donor_mate_ids(vudReadIds)
        self.assertListEqual(self.variantContext.get_unmapped_donor_mate_ids(), vudReadIds, f"The variant context unmapped donor mate ids should have been {vudReadIds}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OTHER GETTER METHODS====================
    def test_getVariantContextLength(self):
        self.assertEqual(self.variantContext.get_variant_context_length(), (abs(self.contextEndAnswer - self.contextStartAnswer)), f"The length of the variant context should have been {abs(self.contextEndAnswer - self.contextStartAnswer)}")
    
    def test_getStartDistanceFromOrigin(self):
        self.assertEqual(self.variantContext.get_start_distance_from_origin(), (abs(self.contextOriginAnswer - self.contextStartAnswer)), f"The start distance from the origin should have been {abs(self.contextOriginAnswer - self.contextStartAnswer)}")
    
    def test_getEndDistanceFromOrigin(self):
        self.assertEqual(self.variantContext.get_end_distance_from_origin(), (abs(self.contextEndAnswer - self.contextOriginAnswer)), f"The end distance from the origin should have been {abs(self.contextEndAnswer - self.contextOriginAnswer)}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT ACCEPTOR READ DATA====================
    def test_getNumberOfVariantContextAcceptorReads(self):
        self.assertEqual(self.variantContext.get_number_of_acceptor_reads(), len(self.acceptorReadsAnswer), f"The variant context acceptor reads are not what was expected")
    
    def test_getVariantContextAcceptorReadIds(self):
        self.assertListEqual(self.variantContext.get_acceptor_read_ids(), self.acceptorReadIdsAnswer, f"The variant context acceptor read ids should have been {self.acceptorReadIdsAnswer}")
    
    def test_getVariantContextAcceptorReadStarts(self):
        self.assertListEqual(self.variantContext.get_acceptor_read_starts(), self.acceptorReadStartsAnswer, f"The variant context acceptor read starts should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorReadLeftPositions(self):
        self.assertListEqual(self.variantContext.get_acceptor_read_left_positions(), self.acceptorReadStartsAnswer, f"The variant context acceptor read left most positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorReadEnds(self):
        self.assertListEqual(self.variantContext.get_acceptor_read_ends(), self.acceptorReadEndsAnswer, f"The variant context acceptor read end positions should have been {self.acceptorReadEndsAnswer}")
    
    def test_getAcceptorReadRightPositions(self):
        self.assertListEqual(self.variantContext.get_acceptor_read_right_positions(), self.acceptorRightPosAnswer, f"The variant context acceptor read right most positions should have been {self.acceptorRightPosAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT DONOR READ DATA====================
    def test_getNumberOfVariantContextDonorReads(self):
        self.assertEqual(self.variantContext.get_number_of_donor_reads(), len(self.donorReadsAnswer), "The variant context donor reads are not what was expected")
    
    def test_getVariantContextDonorReadIds(self):
        self.assertListEqual(self.variantContext.get_donor_read_ids(), self.donorReadIdsAnswer, f"the variant context donor read ids should have been {self.donorReadIdsAnswer}")
    
    def test_getVariantContextDonorReadStarts(self):
        self.assertListEqual(self.variantContext.get_donor_read_starts(), self.donorReadStartsAnswer, f"The variant context donor read start positions should have been {self.donorReadStartsAnswer}")
    
    def test_getVariantContextDonorReadLeftPositions(self):
        self.assertListEqual(self.variantContext.get_donor_read_left_positions(), self.donorLeftPosAnswer, f"The variant context donor read left most positions should have been {self.donorLeftPosAnswer}")
    
    def test_getVariantContextDonorReadEnds(self):
        self.assertListEqual(self.variantContext.get_donor_read_ends(), self.donorReadEndsAnswer, f"The variant context donor read end positions should have been {self.donorReadEndsAnswer}")
    
    def test_getVariantContextDonorReadRightPositions(self):
        self.assertListEqual(self.variantContext.get_donor_read_right_positions(), self.donorReadEndsAnswer, f"The variant context donor read right most positions should have been {self.donorReadEndsAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR ADDING DATA TO THE VARIANT CONTEXT====================
    def test_setAcceptorContext(self):
        accConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.variantContext.set_acceptor_context(accConToAdd)
        self.assertEqual(self.variantContext.get_acceptor_context().to_string(), accConToAdd.to_string(), "The saved acceptor context is not the same as the one which was just set")
    
    def test_setDonorContext(self):
        donConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.variantContext.set_donor_context(donConToAdd)
        self.assertEqual(self.variantContext.get_donor_context().to_string(), donConToAdd.to_string(), "The saved donor context is not the same as the one which was just set")
    
    def test_addAcceptorContext(self):
        accConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.variantContext.add_acceptor_context('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.assertEqual(self.variantContext.get_acceptor_context().to_string(), accConToAdd.to_string(), "The saved donor context is not the same as the one which was just set")
    
    def test_addDonorContext(self):
        donConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.variantContext.add_donor_context('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.assertEqual(self.variantContext.get_donor_context().to_string(), donConToAdd.to_string(), "The saved donor context is not the same as the one which was just set")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING ACCEPTOR CONTEXT DATA====================
    def test_getAcceptorContextId(self):
        self.assertEqual(self.variantContext.get_acceptor_context_id(), self.contextIdAnswer, f"The acceptor context id should have been")
    
    def test_getAcceptorSampleId(self):
        self.assertEqual(self.variantContext.get_acceptor_context_sample_id(), self.contextSampleAnswer, f"The acceptor context sample id should have been {self.contextSampleAnswer}")
    
    def test_getAcceptorContextChrom(self):
        self.assertEqual(self.variantContext.get_acceptor_context_chrom(), self.contextChromAnswer, f"The acceptor context chromosome should have been {self.contextChromAnswer}")
    
    def test_getAcceptorContextOrigin(self):
        self.assertEqual(self.variantContext.get_acceptor_context_origin(), self.contextOriginAnswer, f"The acceptor context origin should have been {self.contextOriginAnswer}")
    
    def test_getAcceptorContextStart(self):
        self.assertEqual(self.variantContext.get_acceptor_context_start(), self.acceptorStartAnswer, f"The acceptor context start positions should have been {self.acceptorStartAnswer}")
    
    def test_getAcceptorContextEnd(self):
        self.assertEqual(self.variantContext.get_acceptor_context_end(), self.acceptorEndAnswer, f"The acceptor context end positions should have been {self.acceptorEndAnswer}")
    
    def test_getAcceptorContextLength(self):
        self.assertEqual(self.variantContext.get_acceptor_context_length(), abs(self.acceptorEndAnswer - self.acceptorStartAnswer), f"The acceptor context length should have been {abs(self.acceptorEndAnswer - self.acceptorStartAnswer)}")
    
    def test_getAcceptorContextReads(self):
        self.assertEqual(self.variantContext.get_acceptor_context_reads(), self.acceptorReadsAnswer, f"The acceptor context reads are not what was expected")
    
    def test_getAcceptorContextReadIds(self):
        self.assertListEqual(self.variantContext.get_acceptor_context_read_ids(), self.acceptorReadIdsAnswer, f"The acceptor context read ids should have been {self.acceptorReadIdsAnswer}")
    
    def test_getAcceptorContextReadStarts(self):
        self.assertEqual(self.variantContext.get_acceptor_context_read_starts(), self.acceptorReadStartsAnswer, f"The acceptor context read start positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorContextReadLeftPositions(self):
        self.assertEqual(self.variantContext.get_acceptor_context_read_left_positions(), self.acceptorReadStartsAnswer, f"The acceptor context read left most positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorContextReadEnds(self):
        self.assertEqual(self.variantContext.get_acceptor_context_read_ends(), self.acceptorReadEndsAnswer, f"The acceptor read end positions should have been {self.acceptorReadEndsAnswer}")
    
    def test_getAcceptorContextReadRightPositions(self):
        self.assertEqual(self.variantContext.get_acceptor_context_read_right_positions(), self.acceptorRightPosAnswer, f"The acceptor context read right most positions should have been {self.acceptorRightPosAnswer}")
    
    def test_getAcceptorContextReadLengths(self):
        self.assertEqual(self.variantContext.get_acceptor_context_read_lengths(), self.readLensAnswer, f"The acceptor context read lengths should have been {self.readLensAnswer}")
    
    def test_getAcceptContextUnmappedMateIds(self):
        auReadIds = ['auRead1', 'auRead2', 'auRead3']
        self.variantContext.set_acceptor_context_unmapped_mates(auReadIds)
        self.assertListEqual(self.variantContext.get_acceptor_context_unmapped_mate_ids(), auReadIds, f"The acceptor context unmapped ids should have been {auReadIds}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING DONOR CONTEXT DATA====================
    def test_getDonorContextId(self):
        self.assertEqual(self.variantContext.get_donor_context_id(), self.contextIdAnswer, f"The donor context id should have been {self.contextIdAnswer}")
    
    def test_getDonorSampleId(self):
        self.assertEqual(self.variantContext.get_donor_context_sample_id(), self.contextSampleAnswer, f"The donor context sample id should have been {self.contextSampleAnswer}")
    
    def test_getDonorContextChrom(self):
        self.assertEqual(self.variantContext.get_donor_context_chrom(), self.contextChromAnswer, f"The donor context chromosome should have been {self.contextChromAnswer}")
    
    def test_getDonorContextOrigin(self):
        self.assertEqual(self.variantContext.get_donor_context_origin(), self.contextOriginAnswer, f"The donor context origin should have been {self.contextOriginAnswer}")
    
    def test_getDonorContextStart(self):
        self.assertEqual(self.variantContext.get_donor_context_start(), self.donorStartAnswer, f"The donor context start position should have been {self.donorStartAnswer}")
    
    def test_getDonorContextEnd(self):
        self.assertEqual(self.variantContext.get_donor_context_end(), self.donorEndAnswer, f"The donor context end position should have been {self.donorEndAnswer}")
    
    def test_getDonorContextLength(self):
        self.assertEqual(self.variantContext.get_donor_context_length(), abs(self.donorEndAnswer - self.donorStartAnswer), f"The donor context length should have been {abs(self.donorEndAnswer - self.donorStartAnswer)}")
    
    def test_getDonorContextReads(self):
        self.assertEqual(self.variantContext.get_donor_context_reads(), self.donorReadsAnswer, f"The donor context reads are not what was expected")
    
    def test_getDonorContextReadIds(self):
        self.assertListEqual(self.variantContext.get_donor_context_read_ids(), self.donorReadIdsAnswer, f"The donor context read ids should have been {self.donorReadIdsAnswer}")
    
    def test_getDonorContextReadStarts(self):
        self.assertEqual(self.variantContext.get_donor_context_read_starts(), self.donorReadStartsAnswer, f"The donor context read start positions should have been {self.donorReadStartsAnswer}")
    
    def test_getDonorContextReadLeftPositions(self):
        self.assertEqual(self.variantContext.get_donor_context_read_left_positions(), self.donorLeftPosAnswer, f"The donor context read left most positions should have been {self.donorLeftPosAnswer}")
    
    def test_getDonorContextReadEnds(self):
        self.assertEqual(self.variantContext.get_donor_context_read_ends(), self.donorReadEndsAnswer, f"The donor context read end positions should have been {self.donorReadEndsAnswer}")
    
    def test_getDonorContextReadRightPositions(self):
        self.assertEqual(self.variantContext.get_donor_context_read_right_positions(), self.donorReadEndsAnswer, f"The donor context read right most positions should have been {self.donorReadEndsAnswer}")
    
    def test_getDonorContextReadLengths(self):
        self.assertEqual(self.variantContext.get_donor_context_read_lengths(), self.donorReadLensAnswer, f"The donor context read lengths should have been {self.readLensAnswer}")
    
    def test_getAcceptContextUnmappedMateIds(self):
        duReadIds = ['duRead1', 'duRead2', 'duRead3']
        self.variantContext.set_donor_context_unmapped_mates(duReadIds)
        self.assertListEqual(self.variantContext.get_donor_context_unmapped_mate_ids(), duReadIds, f"The donor context id should have been {duReadIds}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OUTPUT METHODS====================
    def test_toString(self):
        toStringAnswer = f"{self.contextIdAnswer}\t{self.contextSampleAnswer}\t{self.contextChromAnswer}\t{self.contextOriginAnswer}\t{self.contextStartAnswer}\t{self.contextEndAnswer}\t{abs(self.acceptorEndAnswer - self.acceptorStartAnswer)}\t{abs(self.donorEndAnswer - self.donorStartAnswer)}\t{len(self.acceptorReadsAnswer)}\t{len(self.donorReadsAnswer)}\t{float(len(self.acceptorReadsAnswer)/len(self.donorReadsAnswer))}\t" +';'.join(self.acceptorReadIdsAnswer)+ "\t" +';'.join(self.donorReadIdsAnswer)
        self.assertEqual(self.variantContext.to_string(), toStringAnswer, f"The toString line should have been {toStringAnswer}")
    
    def test_toStatisticsString(self):
        toStatisticsAnswer = f"{self.contextIdAnswer}\t"
        self.assertEqual(self.variantContext.to_statistics_string(), self.toStatisticsAnswer, f"The statistics line should have been {self.toStatisticsAnswer}")
