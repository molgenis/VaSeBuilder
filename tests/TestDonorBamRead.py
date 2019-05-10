import unittest
from DonorBamRead import DonorBamRead

class TestDonorBamRead(unittest.TestCase):
	# Creates the answer variables and the DonorBamRead used for testing
	def setUp(self):
		self.readIdAnswer = 'HHKY2CCXX160108:1:2122:24160:2522'
		self.readPnAnswer = '1'
		self.readChromAnswer = '21'
		self.readPosAnswer = 9411193
		self.readLenAnswer = 151
		self.readSeqAnswer = 'AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGACTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA'
		self.readQualsAnswer = '><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>>==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>'
		self.readMapQAnswer = 40
		self.dbamRead = DonorBamRead(self.readIdAnswer, self.readPnAnswer, self.readChromAnswer, self.readPosAnswer, self.readLenAnswer, self.readSeqAnswer, self.readQualsAnswer, self.readMapQAnswer)
		
		# Create the other variables saving the correct answers that should be returned (makes it easy to modify the test)
		self.qscoresAnswer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
		self.readEndAnswer = 9411344
		self.avgQScoreAnswer = 28.490066225165563
		self.medQScoreAnswer = 28
		self.toStringAnswer = self.readIdAnswer +"\t"+ self.readPnAnswer +"\t"+ self.readChromAnswer +"\t"+ str(self.readPosAnswer) +"\t"+ str(self.readLenAnswer) +"\t"+ self.readSeqAnswer +"\t"+ self.readQualsAnswer +"\t"+ str(self.readMapQAnswer)
		self.fastqPnAnswer = self.readIdAnswer +"/"+ self.readPnAnswer +"\n"+ self.readSeqAnswer +"\n+\n"+ self.readQualsAnswer +"\n"
		self.fastqAnswer = self.readIdAnswer +"\n"+ self.readSeqAnswer +"\n+\n"+ self.readQualsAnswer +"\n"
	
	
	
	# ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
	def test_getBamReadId(self):
		self.assertEqual(self.dbamRead.getBamReadId(), self.readIdAnswer, f"Both BAM read identifiers should have been {self.readIdAnswer}")
	
	def test_getBamReadPairNumber(self):
		self.assertEquals(self.dbamRead.getBamReadPairNumber(), self.readPnAnswer, f"Both read pair numbers should have been {self.readPnAnswer}")
	
	def test_getBamReadChrom(self):
		self.assertEquals(self.dbamRead.getBamReadChrom(), self.readChromAnswer, f"Both read chromosomes should have been {self.readChromAnswer}")
	
	def test_getBamReadRefPos(self):
		self.assertEqual(self.dbamRead.getBamReadRefPos(), self.readPosAnswer, f"Both read positions should have been {self.readPosAnswer}")
	
	def test_getBamReadLength(self):
		self.assertEqual(self.dbamRead.getBamReadLength(), self.readLenAnswer, f"Both read lengths should have been {self.readLenAnswer}")
	
	def test_getBamReadRefEnd(self):
		self.assertEqual(self.dbamRead.getBamReadRefEnd(), self.readEndAnswer, f"Both read end positions should have been {self.readEndAnswer}")
	
	def test_getBamReadQual(self):
		self.assertEqual(self.dbamRead.getBamReadQual(), self.readQualsAnswer, f"Both read qualities should have been {self.readQualsAnswer}")
	
	def test_getBamReadQScores(self):
		self.assertEqual(self.dbamRead.getBamReadQScores(), self.qscoresAnswer, f"Both Q-scores should have been {self.qscoresAnswer}")
	
	def test_getBamReadMapQ(self):
		self.assertEqual(self.dbamRead.getMappingQual(), self.readMapQAnswer, f"Both MapQ values should have been {self.readMapQAnswer}")
	
	
	
	# ====================PERFORM THE TESTS FOR THE STATISTICS METHODS====================
	def test_getAverageQScore(self):
		self.assertEqual(self.dbamRead.getAverageQscore(), self.avgQScoreAnswer, f"Both average Q-scores should have been {self.avgQScoreAnswer}")
	
	def test_getMedianQScore(self):
		self.assertEqual(self.dbamRead.getMedianQScore(), self.medQScoreAnswer, f"Both median Q-scores should have been {self.medQScoreAnswer}")
	
	
	
	# ====================PERFORM THE TESTS FOR CHECKING WHETHER THE READ IS R1 OR R2====================
	def test_isRead1(self):
		self.assertTrue(self.dbamRead.isRead1(), "Donor BAM read should have been read 1")
	
	def test_isRead2(self):
		self.assertFalse(self.dbamRead.isRead2(), "Donor BAM read should not have been read 2")
	
	
	
	# ====================PERFORM THE TESTS FOR RETURNING BAM READ STRING REPRESENTATIONS====================
	def test_toString(self):
		self.assertEqual
	
	def test_getAsFastQSeq_pairnum(self):
		self.assertEqual(self.dbamRead.getAsFastQSeq(True), "@"+self.fastqPnAnswer, f"Both answers should have been {self.fastqPnAnswer}")
	
	def test_getAsFastQSeq_nopairnum(self):
		self.assertEqual(self.dbamRead.getAsFastQSeq(), "@"+self.fastqAnswer, f"Both answers should have been {self.fastqAnswer}")
