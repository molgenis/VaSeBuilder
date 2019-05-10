import unittest
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext

class TestOverlapContext(unittest.TestCase):
	# Creates the required 
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
		self.contextStartDistance = self.contextOriginAnswer - self.contextStartAnswer
		self.contextEndDistance = self.contextEndAnswer - self.contextOriginAnswer
		self.overlapContext = OverlapContext(self.contextIdAnswer, self.contextSampleAnswer, self.contextChromAnswer, self.contextOriginAnswer, self.contextStartAnswer, self.contextEndAnswer, self.contextBamReadsAnswer)
		
		# Create the variables containing the context read test answers
		self.numOfContextReads = len(self.contextBamReadsAnswer)
		self.contextReadIdsAnswer = ['HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522', 'HHKY2CCXX160108:1:2122:24160:2522']
		self.contextStartsAnswer = [9411193, 9411193, 9411193]
		self.contextLeftPosAnswer = [9411193, 9411193, 9411193]
		self.contextEndsAnswer = [9411344, 9411344, 9411344]
		self.contextReadLensAnswer = [151,151,151]
		self.contextRightPosAnswer = []
		self.contextReadSeqsAnswer = [self.readSeqAnswer, self.readSeqAnswer, self.readSeqAnswer]
		self.contextReadQualsAnswer = [self.readQualsAnswer, self.readQualsAnswer, self.readQualsAnswer]
		self.qscoreAnswer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
		self.contextReadQScoresAnswer = self.qscoreAnswer + self.qscoreAnswer + self.qscoreAnswer
		self.contextReadMapQsAnswer = [40, 40, 40]
		self.readInContextAnswer_neg = 'HHKY2CCXX160108:1:2122:24160:2555'
		
		# Create the variables containing the context statistics answers.
		self.avgMedReadLenAnswer = [151, 151]
		self.avgMedReadQualAnswer = [28.490066225165563, 28.490066225165563]
		self.avgMedReadMapQAnswer = [40, 40]
		
		# Create the variables containing the overlap context other methods answers
		self.toStringAnswer = f"{self.contextIdAnswer}\t{self.contextSampleAnswer}\t{self.contextChromAnswer}\t{self.contextOriginAnswer}\t{self.contextStartAnswer}\t{self.contextEndAnswer}\t{len(self.contextBamReadsAnswer)}\t" +';'.join([x.getBamReadId() for x in self.contextBamReadsAnswer])
		self.toStatisticsAnswer = f"{self.contextIdAnswer}\t{self.avgMedReadLenAnswer[0]}\t{self.avgMedReadLenAnswer[1]}\t{self.avgMedReadQualAnswer[0]}\t{self.avgMedReadQualAnswer[1]}\t{self.avgMedReadMapQAnswer[0]}\t{self.avgMedReadMapQAnswer[1]}"
		self.compareAnswer = {}
	
	
	
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
		self.assertEqual(self.overlapContext.getStartDistanceFromOrigin(), self.contextStartDistance, f"The start distance from the context origin should have been {self.contextStartDistance}")
	
	def test_getEndDistanceFromOrigin(self):
		self.assertEqual(self.overlapContext.getEndDistanceFromOrigin(), self.contextEndDistance, f"The end distance from the context origin should have been {self.contextEndDistance}")
	
	
	
	# ====================PERFORM THE TESTS FOR GETTING CONTEXT READ INFO====================
	def test_getNumberOfContextReads(self):
		self.assertEqual(self.overlapContext.getNumberOfContextReads(), self.numOfContextReads, f"The number of context reads should have been {self.numOfContextReads}")
	
	def test_getContextBamReadIds(self):
		self.assertEqual(self.overlapContext.getContextBamReadIds(), self.contextReadIdsAnswer, f"The list of context read ids should have been {self.contextReadIdsAnswer}")
	
	def test_getContextBamReadStarts(self):
		self.assertListEqual(self.overlapContext.getContextBamReadStarts(), self.contextStartsAnswer, f"The list of start positions should have been {self.contextStartsAnswer}")
	
	def test_getContextBamReadLeftPositions(self):
		self.assertEqual(self.overlapContext.getContextBamReadLeftPositions(), self.contextLeftPosAnswer, f"The list of left most positions should have been {self.contextLeftPosAnswer}")
	
	def test_getContextBamReadEnds(self):
		self.assertListEqual(self.overlapContext.getContextBamReadEnds(), self.contextEndsAnswer, f"The list of end positions should have been {self.contextEndsAnswer}")
	
	def test_getContextBamReadRightPositions(self):
		self.assertEqual(self.overlapContext.getContextBamReadRightPositions(), self.contextRightPosAnswer, f"The list of right most positions should have been {self.contextRightPosAnswer}")
	
	def test_getContextBamReadLength(self):
		self.assertEqual(self.overlapContext.getContextBamReadLength(), self.contextReadLensAnswer, f"The context read lengths should have been {self.contextReadLensAnswer}")
	
	def test_getContextBamReadSeqs(self):
		self.assertListEqual(self.overlapContext.getContextBamReadSeqs(), self.contextReadSeqsAnswer, f"The context read sequences should have been {self.contextReadSeqsAnswer}")
	
	def test_getContextBamReadQualities(self):
		self.assertEqual(self.overlapContext.getContextBamReadQualities(), self.contextReadQualsAnswer, f"The read qualites should have been {self.contextReadQualsAnswer}")
	
	#def test_getContextBamReadQScores(self):
		#self.assertListEqual(self.overlapContext.getContextBamReadQScores(), self.contextReadQScoresAnswer, f"The Q-scores should have been {self.contextReadQScoresAnswer}")
	
	def test_getContextBamReadMapQs(self):
		self.assertListEqual(self.overlapContext.getContextBamReadMapQs(), self.contextReadMapQsAnswer, f"The qualities should have been {self.contextReadMapQsAnswer}")
	
	def test_readIsInContext_pos(self):
		self.assertTrue(self.overlapContext.readIsInContext(self.readIdAnswer), f"Read id {self.readIdAnswer} should have been in the context")
	
	def test_readIsInContext_neg(self):
		self.assertFalse(self.overlapContext.readIsInContext(self.readInContextAnswer_neg), f"Read id {self.readInContextAnswer_neg} should not have been in the context")
	
	
	
	# ====================PERFORM THE TESTS FOR THE OVERLAP CONTEXT STATISTICS====================
	def test_getAverageAndMedianReadLength(self):
		self.assertListEqual(self.overlapContext.getAverageAndMedianReadLength(), self.avgMedReadLenAnswer, f"Average and median read length should have been {self.avgMedReadLenAnswer}")
	
	def test_getAverageAndMedianReadQual(self):
		self.assertListEqual(self.overlapContext.getAverageAndMedianReadQual(), self.avgMedReadQualAnswer, f"Average and median read length should have been {self.avgMedReadQualAnswer}")
	
	def test_getAverageAndMedianReadMapQ(self):
		self.assertListEqual(self.overlapContext.getAverageAndMedianReadMapQ(), self.avgMedReadMapQAnswer, f"Average and median read length should have been {self.avgMedReadMapQAnswer}")
	
	
	
	# ====================PERFOM THE TESTS FOR THE OTHER METHODS====================
	def test_toString(self):
		self.assertEqual(self.overlapContext.toString(), self.toStringAnswer, f"String representation should have been {self.toStringAnswer}")
	
	def test_toStatisticsString(self):
		self.assertEqual(self.overlapContext.toStatisticsString(), self.toStatisticsAnswer, f"Statistics string should have been {self.toStatisticsAnswer}")
	
	def test_compare(self):
		self.assertDictEqual(self.overlapContext.compare(self.overlapContext), self.compareAnswer, f"Compare results dictionary should have been {self.compareAnswer}")
		
