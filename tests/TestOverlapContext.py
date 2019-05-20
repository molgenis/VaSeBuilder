import unittest
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext

class TestOverlapContext(unittest.TestCase):
    # Creates the variables that are needed for each test method
    def setUp(self):
        # Create three DonorBamRead objects to add to the overlap context
        self.readIdAnswer = 'HHKY2CCXX160108:1:2122:24160:2522'
        self.readPnAnswer = '1'
        self.readChromAnswer = '21'
        self.readPosAnswer = 9411193
        self.readLenAnswer = 151
        self.readSeqAnswer = 'AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGACTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA'
        self.readQualsAnswer = '><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>>==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>'
        self.readMapQAnswer = 40
        self.conRead1 = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
        self.conRead2 = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
        self.conRead3 = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
        
        # Create the overlap context to test
        self.contextIdAnswer = '21_9411259'
        self.contextSampleAnswer = 'testSample'
        self.contextChromAnswer = '21'
        self.contextOriginAnswer = 9411259
        self.contextStartAnswer = 9411193
        self.contextEndAnswer = 9411344
        self.contextLenAnswer = 151
        self.contextBamReadsAnswer = [self.conRead1, self.conRead2, self.conRead3]
        self.unmappedAnswer = []
        self.overlapContext = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.contextStartAnswer, self.contextEndAnswer, self.contextBamReadsAnswer)
        
        # Create the variables containing the context read test answers
        self.qscoreAnswer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
        self.avgMedReadLenAnswer = [151, 151]
        self.avgMedReadQualAnswer = [28.490066225165563, 28.490066225165563]
        self.avgMedReadMapQAnswer = [40, 40]
    
    
    
    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_getContextId(self):
        self.assertEqual(self.overlapContext.getContextId(), self.contextIdAnswer, f"The context id should have been {self.contextIdAnswer}")
    
    def test_getSampleId(self):
        self.assertEqual(self.overlapContext.getSampleId(), self.contextSampleAnswer, f"The sample id should have been {self.contextSampleAnswer}")
    
    def test_getContextChrom(self):
        self.assertEqual(self.overlapContext.getContextChrom(), self.contextChromAnswer, f"The context chromosome should have been {self.contextChromAnswer}")
    
    def test_getContextOrigin(self):
        self.assertEqual(self.overlapContext.getContextOrigin(), self.contextOriginAnswer, f"The context origin should have been {self.contextOriginAnswer}")
    
    def test_getContextStart(self):
        self.assertEqual(self.overlapContext.getContextStart(), self.contextStartAnswer, f"The context start position should have been {self.contextStartAnswer}")
    
    def test_getContextEnd(self):
        self.assertEqual(self.overlapContext.getContextEnd(), self.contextEndAnswer, f"Both context end positions should have been {self.contextEndAnswer}")
    
    def test_getContextBamReads(self):
        self.assertListEqual(self.overlapContext.getContextBamReads(), self.contextBamReadsAnswer, f"The BAM read list should have been {self.contextBamReadsAnswer}")
    
    def test_getUnmappedReadMateIds(self):
        self.assertListEqual(self.overlapContext.getUnmappedReadMateIds(), self.unmappedAnswer, f"The unmapped mate read id list should have been {self.unmappedAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OTHER GETTER METHODS====================
    def test_getContextLength(self):
        self.assertEqual(self.overlapContext.getContextLength(), self.contextLenAnswer, f"The context length should have been {self.contextLenAnswer}")
    
    def test_getStartDistanceFromOrigin(self):
        contextStartDistanceAnswer = self.contextOriginAnswer - self.contextStartAnswer
        self.assertEqual(self.overlapContext.getStartDistanceFromOrigin(), contextStartDistanceAnswer, f"The start distance from the context origin should have been {contextStartDistanceAnswer}")
    
    def test_getEndDistanceFromOrigin(self):
        contextEndDistanceAnswer = self.contextEndAnswer - self.contextOriginAnswer
        self.assertEqual(self.overlapContext.getEndDistanceFromOrigin(), contextEndDistanceAnswer, f"The end distance from the context origin should have been {contextEndDistanceAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR GETTING CONTEXT READ INFO====================
    def test_getNumberOfContextReads(self):
        numOfContextReadsAnswer = len(self.contextBamReadsAnswer)
        self.assertEqual(self.overlapContext.getNumberOfContextReads(), numOfContextReadsAnswer, f"The number of context reads should have been {numOfContextReadsAnswer}")
    
    def test_getContextBamReadIds(self):
        contextReadIdsAnswer = ['HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522']
        self.assertEqual(self.overlapContext.getContextBamReadIds(), contextReadIdsAnswer, f"The list of context read ids should have been {contextReadIdsAnswer}")
    
    def test_getContextBamReadStarts(self):
        contextStartsAnswer = [9411193, 9411193, 9411193]
        self.assertListEqual(self.overlapContext.getContextBamReadStarts(), contextStartsAnswer, f"The list of start positions should have been {contextStartsAnswer}")
    
    def test_getContextBamReadLeftPositions(self):
        contextLeftPosAnswer = [9411193, 9411193, 9411193]
        self.assertEqual(self.overlapContext.getContextBamReadLeftPositions(), contextLeftPosAnswer, f"The list of left most positions should have been {contextLeftPosAnswer}")
    
    def test_getContextBamReadEnds(self):
        contextEndsAnswer = [9411344, 9411344, 9411344]
        self.assertListEqual(self.overlapContext.getContextBamReadEnds(), contextEndsAnswer, f"The list of end positions should have been {contextEndsAnswer}")
    
    def test_getContextBamReadRightPositions(self):
        contextRightPosAnswer = []
        self.assertEqual(self.overlapContext.getContextBamReadRightPositions(), contextRightPosAnswer, f"The list of right most positions should have been {contextRightPosAnswer}")
    
    def test_getContextBamReadLengths(self):
        contextReadLensAnswer = [151,151,151]
        self.assertEqual(self.overlapContext.getContextBamReadLengths(), contextReadLensAnswer, f"The context read lengths should have been {contextReadLensAnswer}")
    
    def test_getContextBamReadSeqs(self):
        contextReadSeqsAnswer = [self.readSeqAnswer, self.readSeqAnswer, self.readSeqAnswer]
        self.assertListEqual(self.overlapContext.getContextBamReadSeqs(), contextReadSeqsAnswer, f"The context read sequences should have been {contextReadSeqsAnswer}")
    
    def test_getContextBamReadQualities(self):
        contextReadQualsAnswer = [self.readQualsAnswer, self.readQualsAnswer, self.readQualsAnswer]
        self.assertEqual(self.overlapContext.getContextBamReadQualities(), contextReadQualsAnswer, f"The read qualites should have been {contextReadQualsAnswer}")
    
    #def test_getContextBamReadQScores(self):
        #contextReadQScoresAnswer = self.qscoreAnswer + self.qscoreAnswer + self.qscoreAnswer
        #self.assertListEqual(self.overlapContext.getContextBamReadQScores(), contextReadQScoresAnswer, f"The Q-scores should have been {contextReadQScoresAnswer}")
    
    def test_getContextBamReadMapQs(self):
        contextReadMapQsAnswer = [40, 40, 40]
        self.assertListEqual(self.overlapContext.getContextBamReadMapQs(), contextReadMapQsAnswer, f"The qualities should have been {contextReadMapQsAnswer}")
    
    def test_readIsInContext_pos(self):
        self.assertTrue(self.overlapContext.readIsInContext(self.readIdAnswer), f"Read id {self.readIdAnswer} should have been in the context")
    
    def test_readIsInContext_neg(self):
        readInContextAnswer_neg = 'HHKY2CCXX160108:1:2122:24160:2555'
        self.assertFalse(self.overlapContext.readIsInContext(readInContextAnswer_neg), f"Read id {readInContextAnswer_neg} should not have been in the context")
    
    
    
    # ====================PERFORM THE TESTS FOR THE OVERLAP CONTEXT STATISTICS====================
    def test_getAverageAndMedianReadLength(self):
        self.assertListEqual(self.overlapContext.getAverageAndMedianReadLength(), self.avgMedReadLenAnswer, f"Average and median read length should have been {self.avgMedReadLenAnswer}")
    
    def test_getAverageAndMedianReadQual(self):
        self.assertListEqual(self.overlapContext.getAverageAndMedianReadQual(), self.avgMedReadQualAnswer, f"Average and median read length should have been {self.avgMedReadQualAnswer}")
    
    def test_getAverageAndMedianReadMapQ(self):
        self.assertListEqual(self.overlapContext.getAverageAndMedianReadMapQ(), self.avgMedReadMapQAnswer, f"Average and median read length should have been {self.avgMedReadMapQAnswer}")
    
    
    
    # ====================PERFOM THE TESTS FOR THE OTHER METHODS====================
    def test_toString(self):
        toStringAnswer = f"{self.contextIdAnswer}\t{self.contextSampleAnswer}\t{self.contextChromAnswer}\t{self.contextOriginAnswer}\t{self.contextStartAnswer}\t{self.contextEndAnswer}\t{len(self.contextBamReadsAnswer)}\t" +';'.join([x.getBamReadId() for x in self.contextBamReadsAnswer])
        self.assertEqual(self.overlapContext.toString(), toStringAnswer, f"String representation should have been {toStringAnswer}")
    
    def test_toStatisticsString(self):
        toStatisticsAnswer = f"{self.contextIdAnswer}\t{self.avgMedReadLenAnswer[0]}\t{self.avgMedReadLenAnswer[1]}\t{self.avgMedReadQualAnswer[0]}\t{self.avgMedReadQualAnswer[1]}\t{self.avgMedReadMapQAnswer[0]}\t{self.avgMedReadMapQAnswer[1]}"
        self.assertEqual(self.overlapContext.toStatisticsString(), toStatisticsAnswer, f"Statistics string should have been {toStatisticsAnswer}")
    
    def test_compare(self):
        compareAnswer = {}
        self.assertDictEqual(self.overlapContext.compare(self.overlapContext), compareAnswer, f"Compare results dictionary should have been {compareAnswer}")
