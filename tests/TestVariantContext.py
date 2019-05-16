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
        self.toStringAnswer = f"{self.contextIdAnswer}\t{self.contextSampleAnswer}\t{self.contextChromAnswer}\t{self.contextOriginAnswer}\t{self.contextStartAnswer}\t{self.contextEndAnswer}\t{abs(self.acceptorEndAnswer - self.acceptorStartAnswer)}\t{abs(self.donorEndAnswer - self.donorStartAnswer)}\t{len(self.acceptorReadsAnswer)}\t{len(self.donorReadsAnswer)}\t{float(len(self.acceptorReadsAnswer)/len(self.donorReadsAnswer))}\t" +';'.join(self.acceptorReadIdsAnswer)+ "\t" +';'.join(self.donorReadIdsAnswer)
        self.toStatisticsAnswer = f"{self.contextIdAnswer}\t"
    
    
    
    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_getVariantContextId(self):
        self.assertEqual(self.variantContext.getVariantContextId(), self.contextIdAnswer, f"The variant context id should have been {self.contextIdAnswer}")
    
    def test_getVariantContextSample(self):
        self.assertEqual(self.variantContext.getVariantContextSample(), self.contextSampleAnswer, f"The variant context sample should have been {self.contextSampleAnswer}")
    
    def test_getVariantContextChrom(self):
        self.assertEqual(self.variantContext.getVariantContextChrom(), self.contextChromAnswer, f"The variant context chromosome should have been {self.contextChromAnswer}")
    
    def test_getVariantContextOrigin(self):
        self.assertEqual(self.variantContext.getVariantContextOrigin(), self.contextOriginAnswer, f"The variant context origin position should have been {self.contextOriginAnswer}")
    
    def test_getVariantContextStart(self):
        self.assertEqual(self.variantContext.getVariantContextStart(), self.contextStartAnswer, f"The variant context start position should have been {self.contextStartAnswer}")
    
    def test_getVariantContextEnd(self):
        self.assertEqual(self.variantContext.getVariantContextEnd(), self.contextEndAnswer, f"The variant context end position should have been {self.contextEndAnswer}")
    
    def test_getVariantContextAcceptorReads(self):
        self.assertListEqual(self.variantContext.getAcceptorReads(), self.acceptorReadsAnswer, f"The variant context acceptor reads are not what was expected")
    
    def test_getVariantContextDonorReads(self):
        self.assertListEqual(self.variantContext.getDonorReads(), self.donorReadsAnswer, f"The variant context donor reads are not what was expected")
    
    def test_getAcceptorContext(self):
        self.assertEqual(self.variantContext.getAcceptorContext(), self.acceptorContextAnswer, f"The acceptor context is not what was expected")
    
    def test_getDonorContext(self):
        self.assertEqual(self.variantContext.getDonorContext(), self.donorContextAnswer, f"The donor context is not what was expected")
    
    #def test_getUnmappedAcceptorMateIds(self):
        #self.assertListEqual(self.variantContext.getUnmappedAcceptorMateIds(), self., f"The variant context unmapped acceptor mate ids should have been {}")
    
    #def test_getUnmappedDonorMateIds(self):
        #self.assertListEqual(self.variantContext.getUnmappedDonorMateIds(), self., f"The variant context unmapped donor mate ids should have been {}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OTHER GETTER METHODS====================
    def test_getVariantContextLength(self):
        self.assertEqual(self.variantContext.getVariantContextLength(), (abs(self.contextEndAnswer - self.contextStartAnswer)), f"The length of the variant context should have been {abs(self.contextEndAnswer - self.contextStartAnswer)}")
    
    def test_getStartDistanceFromOrigin(self):
        self.assertEqual(self.variantContext.getStartDistanceFromOrigin(), (abs(self.contextOriginAnswer - self.contextStartAnswer)), f"The start distance from the origin should have been {abs(self.contextOriginAnswer - self.contextStartAnswer)}")
    
    def test_getEndDistanceFromOrigin(self):
        self.assertEqual(self.variantContext.getEndDistanceFromOrigin(), (abs(self.contextEndAnswer - self.contextOriginAnswer)), f"The end distance from the origin should have been {abs(self.contextEndAnswer - self.contextOriginAnswer)}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT ACCEPTOR READ DATA====================
    def test_getNumberOfVariantContextAcceptorReads(self):
        self.assertEqual(self.variantContext.getNumberOfAcceptorReads(), len(self.acceptorReadsAnswer), f"The variant context acceptor reads are not what was expected")
    
    def test_getVariantContextAcceptorReadIds(self):
        self.assertListEqual(self.variantContext.getAcceptorReadIds(), self.acceptorReadIdsAnswer, f"The variant context acceptor read ids should have been {self.acceptorReadIdsAnswer}")
    
    def test_getVariantContextAcceptorReadStarts(self):
        self.assertListEqual(self.variantContext.getAcceptorReadStarts(), self.acceptorReadStartsAnswer, f"The variant context acceptor read starts should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorReadLeftPositions(self):
        self.assertListEqual(self.variantContext.getAcceptorReadLeftPositions(), self.acceptorReadStartsAnswer, f"The variant context acceptor read left most positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorReadEnds(self):
        self.assertListEqual(self.variantContext.getAcceptorReadEnds(), self.acceptorReadEndsAnswer, f"The variant context acceptor read end positions should have been {self.acceptorReadEndsAnswer}")
    
    def test_getAcceptorReadRightPositions(self):
        self.assertListEqual(self.variantContext.getAcceptorReadRightPositions(), self.acceptorRightPosAnswer, f"The variant context acceptor read right most positions should have been {self.acceptorRightPosAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT DONOR READ DATA====================
    def test_getNumberOfVariantContextDonorReads(self):
        self.assertEqual(self.variantContext.getNumberOfDonorReads(), len(self.donorReadsAnswer), "The variant context donor reads are not what was expected")
    
    def test_getVariantContextDonorReadIds(self):
        self.assertListEqual(self.variantContext.getDonorReadIds(), self.donorReadIdsAnswer, f"the variant context donor read ids should have been {self.donorReadIdsAnswer}")
    
    def test_getVariantContextDonorReadStarts(self):
        self.assertListEqual(self.variantContext.getDonorReadStarts(), self.donorReadStartsAnswer, f"The variant context donor read start positions should have been {self.donorReadStartsAnswer}")
    
    def test_getVariantContextDonorReadLeftPositions(self):
        self.assertListEqual(self.variantContext.getDonorReadLeftPositions(), self.donorLeftPosAnswer, f"The variant context donor read left most positions should have been {self.donorLeftPosAnswer}")
    
    def test_getVariantContextDonorReadEnds(self):
        self.assertListEqual(self.variantContext.getDonorReadEnds(), self.donorReadEndsAnswer, f"The variant context donor read end positions should have been {self.donorReadEndsAnswer}")
    
    def test_getVariantContextDonorReadRightPositions(self):
        self.assertListEqual(self.variantContext.getDonorReadRightPositions(), self.donorReadEndsAnswer, f"The variant context donor read right most positions should have been {self.donorReadEndsAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR ADDING DATA TO THE VARIANT CONTEXT====================
    def test_setAcceptorContext(self):
        accConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.variantContext.setAcceptorContext(accConToAdd)
        self.assertEqual(self.variantContext.getAcceptorContext().toString(), accConToAdd.toString(), "The saved acceptor context is not the same as the one which was just set")
    
    def test_setDonorContext(self):
        donConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.variantContext.setDonorContext(donConToAdd)
        self.assertEqual(self.variantContext.getDonorContext().toString(), donConToAdd.toString(), "The saved donor context is not the same as the one which was just set")
    
    def test_addAcceptorContext(self):
        accConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.variantContext.addAcceptorContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptorReadsAnswer)
        self.assertEqual(self.variantContext.getAcceptorContext().toString(), accConToAdd.toString(), "The saved donor context is not the same as the one which was just set")
    
    def test_addDonorContext(self):
        donConToAdd = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.variantContext.addDonorContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donorReadsAnswer)
        self.assertEqual(self.variantContext.getDonorContext().toString(), donConToAdd.toString(), "The saved donor context is not the same as the one which was just set")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING ACCEPTOR CONTEXT DATA====================
    def test_getAcceptorContextId(self):
        self.assertEqual(self.variantContext.getAcceptorContextId(), self.contextIdAnswer, f"The acceptor context id should have been")
    
    def test_getAcceptorSampleId(self):
        self.assertEqual(self.variantContext.getAcceptorContextSampleId(), self.contextSampleAnswer, f"The acceptor context sample id should have been {self.contextSampleAnswer}")
    
    def test_getAcceptorContextChrom(self):
        self.assertEqual(self.variantContext.getAcceptorContextChrom(), self.contextChromAnswer, f"The acceptor context chromosome should have been {self.contextChromAnswer}")
    
    def test_getAcceptorContextOrigin(self):
        self.assertEqual(self.variantContext.getAcceptorContextOrigin(), self.contextOriginAnswer, f"The acceptor context origin should have been {self.contextOriginAnswer}")
    
    def test_getAcceptorContextStart(self):
        self.assertEqual(self.variantContext.getAcceptorContextStart(), self.acceptorStartAnswer, f"The acceptor context start positions should have been {self.acceptorStartAnswer}")
    
    def test_getAcceptorContextEnd(self):
        self.assertEqual(self.variantContext.getAcceptorContextEnd(), self.acceptorEndAnswer, f"The acceptor context end positions should have been {self.acceptorEndAnswer}")
    
    def test_getAcceptorContextLength(self):
        self.assertEqual(self.variantContext.getAcceptorContextLength(), abs(self.acceptorEndAnswer - self.acceptorStartAnswer), f"The acceptor context length should have been {abs(self.acceptorEndAnswer - self.acceptorStartAnswer)}")
    
    def test_getAcceptorContextReads(self):
        self.assertEqual(self.variantContext.getAcceptorContextReads(), self.acceptorReadsAnswer, f"The acceptor context reads are not what was expected")
    
    def test_getAcceptorContextReadIds(self):
        self.assertListEqual(self.variantContext.getAcceptorContextReadIds(), self.acceptorReadIdsAnswer, f"The acceptor context read ids should have been {self.acceptorReadIdsAnswer}")
    
    def test_getAcceptorContextReadStarts(self):
        self.assertEqual(self.variantContext.getAcceptorContextReadStarts(), self.acceptorReadStartsAnswer, f"The acceptor context read start positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorContextReadLeftPositions(self):
        self.assertEqual(self.variantContext.getAcceptorContextReadLeftPositions(), self.acceptorReadStartsAnswer, f"The acceptor context read left most positions should have been {self.acceptorReadStartsAnswer}")
    
    def test_getAcceptorContextReadEnds(self):
        self.assertEqual(self.variantContext.getAcceptorContextReadEnds(), self.acceptorReadEndsAnswer, f"The acceptor read end positions should have been {self.acceptorReadEndsAnswer}")
    
    def test_getAcceptorContextReadRightPositions(self):
        self.assertEqual(self.variantContext.getAcceptorContextReadRightPositions(), self.acceptorRightPosAnswer, f"The acceptor context read right most positions should have been {self.acceptorRightPosAnswer}")
    
    def test_getAcceptorContextReadLengths(self):
        self.assertEqual(self.variantContext.getAcceptorContextReadLengths(), self.readLensAnswer, f"The acceptor context read lengths should have been {self.readLensAnswer}")
    
    #def test_getAcceptContextUnmappedMateIds(self):
        #self.assertEqual(self.variantContext.(), self., f"The acceptor context id should have been")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING DONOR CONTEXT DATA====================
    def test_getDonorContextId(self):
        self.assertEqual(self.variantContext.getDonorContextId(), self.contextIdAnswer, f"The donor context id should have been {self.contextIdAnswer}")
    
    def test_getDonorSampleId(self):
        self.assertEqual(self.variantContext.getDonorContextSampleId(), self.contextSampleAnswer, f"The donor context sample id should have been {self.contextSampleAnswer}")
    
    def test_getDonorContextChrom(self):
        self.assertEqual(self.variantContext.getDonorContextChrom(), self.contextChromAnswer, f"The donor context chromosome should have been {self.contextChromAnswer}")
    
    def test_getDonorContextOrigin(self):
        self.assertEqual(self.variantContext.getDonorContextOrigin(), self.contextOriginAnswer, f"The donor context origin should have been {self.contextOriginAnswer}")
    
    def test_getDonorContextStart(self):
        self.assertEqual(self.variantContext.getDonorContextStart(), self.donorStartAnswer, f"The donor context start position should have been {self.donorStartAnswer}")
    
    def test_getDonorContextEnd(self):
        self.assertEqual(self.variantContext.getDonorContextEnd(), self.donorEndAnswer, f"The donor context end position should have been {self.donorEndAnswer}")
    
    def test_getDonorContextLength(self):
        self.assertEqual(self.variantContext.getDonorContextLength(), abs(self.donorEndAnswer - self.donorStartAnswer), f"The donor context length should have been {abs(self.donorEndAnswer - self.donorStartAnswer)}")
    
    def test_getDonorContextReads(self):
        self.assertEqual(self.variantContext.getDonorContextReads(), self.donorReadsAnswer, f"The donor context reads are not what was expected")
    
    def test_getDonorContextReadIds(self):
        self.assertListEqual(self.variantContext.getDonorContextReadIds(), self.donorReadIdsAnswer, f"The donor context read ids should have been {self.donorReadIdsAnswer}")
    
    def test_getDonorContextReadStarts(self):
        self.assertEqual(self.variantContext.getDonorContextReadStarts(), self.donorReadStartsAnswer, f"The donor context read start positions should have been {self.donorReadStartsAnswer}")
    
    def test_getDonorContextReadLeftPositions(self):
        self.assertEqual(self.variantContext.getDonorContextReadLeftPositions(), self.donorLeftPosAnswer, f"The donor context read left most positions should have been {self.donorLeftPosAnswer}")
    
    def test_getDonorContextReadEnds(self):
        self.assertEqual(self.variantContext.getDonorContextReadEnds(), self.donorReadEndsAnswer, f"The donor context read end positions should have been {self.donorReadEndsAnswer}")
    
    def test_getDonorContextReadRightPositions(self):
        self.assertEqual(self.variantContext.getDonorContextReadRightPositions(), self.donorReadEndsAnswer, f"The donor context read right most positions should have been {self.donorReadEndsAnswer}")
    
    def test_getDonorContextReadLengths(self):
        self.assertEqual(self.variantContext.getDonorContextReadLengths(), self.donorReadLensAnswer, f"The donor context read lengths should have been {self.readLensAnswer}")
    
    #def test_getAcceptContextUnmappedMateIds(self):
        #self.assertEqual(self.variantContext.(), self., f"The donor context id should have been")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OUTPUT METHODS====================
    def test_toString(self):
        self.assertEqual(self.variantContext.toString(), self.toStringAnswer, f"The toString line should have been {self.toStringAnswer}")
    
    #def test_toStatisticsString(self):
        #self.assertEqual(self.variantContext.toStatisticsString(), self.toStatisticsAnswer, f"The statistics line should have been {self.toStatisticsAnswer}")
