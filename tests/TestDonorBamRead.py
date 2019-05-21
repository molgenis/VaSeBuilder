import unittest
from DonorBamRead import DonorBamRead

class TestDonorBamRead(unittest.TestCase):
    # Creates the answer variables and the DonorBamRead used for testing
    def setUp(self):
        self.read_id_answer = 'HHKY2CCXX160108:1:2122:24160:2522'
        self.read_pn_answer = '1'
        self.read_chrom_answer = '21'
        self.read_pos_answer = 9411193
        self.read_len_answer = 151
        self.read_seq_answer = 'AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGACTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA'
        self.read_quals_answer = '><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>>==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>'
        self.read_map_q_answer = 40
        self.dbam_read = DonorBamRead(self.read_id_answer, self.read_pn_answer, self.read_chrom_answer, self.read_pos_answer, self.read_len_answer, self.read_seq_answer, self.read_quals_answer, self.read_map_q_answer)
    
    
    
    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_getBamReadId(self):
        self.assertEqual(self.dbam_read.get_bam_read_id(), self.read_id_answer, f"Both BAM read identifiers should have been {self.read_id_answer}")
    
    def test_getBamReadPairNumber(self):
        self.assertEquals(self.dbam_read.get_bam_read_pair_number(), self.read_pn_answer, f"Both read pair numbers should have been {self.read_pn_answer}")
    
    def test_getBamReadChrom(self):
        self.assertEquals(self.dbam_read.get_bam_read_chrom(), self.read_chrom_answer, f"Both read chromosomes should have been {self.read_chrom_answer}")
    
    def test_getBamReadRefPos(self):
        self.assertEqual(self.dbam_read.get_bam_read_ref_pos(), self.read_pos_answer, f"Both read positions should have been {self.read_pos_answer}")
    
    def test_getBamReadLength(self):
        self.assertEqual(self.dbam_read.get_bam_read_length(), self.read_len_answer, f"Both read lengths should have been {self.read_len_answer}")
    
    def test_getBamReadRefEnd(self):
        readEndAnswer = 9411344
        self.assertEqual(self.dbam_read.get_bam_read_ref_end(), readEndAnswer, f"Both read end positions should have been {readEndAnswer}")
    
    def test_getBamReadQual(self):
        self.assertEqual(self.dbam_read.get_bam_read_qual(), self.read_quals_answer, f"Both read qualities should have been {self.read_quals_answer}")
    
    def test_getBamReadQScores(self):
        qscoresAnswer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
        self.assertEqual(self.dbam_read.get_bam_read_q_scores(), qscoresAnswer, f"Both Q-scores should have been {qscoresAnswer}")
    
    def test_getBamReadMapQ(self):
        self.assertEqual(self.dbam_read.get_mapping_qual(), self.read_map_q_answer, f"Both MapQ values should have been {self.read_map_q_answer}")
    
    
    
    # ====================PERFORM THE TESTS FOR THE STATISTICS METHODS====================
    def test_getAverageQScore(self):
        avg_qscore_answer = 28.490066225165563
        self.assertEqual(self.dbam_read.get_average_qscore(), avg_qscore_answer, f"Both average Q-scores should have been {avg_qscore_answer}")
    
    def test_getMedianQScore(self):
        medQScoreAnswer = 28
        self.assertEqual(self.dbam_read.get_median_qscore(), medQScoreAnswer, f"Both median Q-scores should have been {medQScoreAnswer}")
    
    
    
    # ====================PERFORM THE TESTS FOR CHECKING WHETHER THE READ IS R1 OR R2====================
    def test_isRead1(self):
        self.assertTrue(self.dbam_read.is_read1(), "Donor BAM read should have been read 1")
    
    def test_isRead2(self):
        self.assertFalse(self.dbam_read.is_read2(), "Donor BAM read should not have been read 2")
    
    
    
    # ====================PERFORM THE TESTS FOR RETURNING BAM READ STRING REPRESENTATIONS====================
    def test_toString(self):
        to_string_answer = f"{self.read_id_answer}\t{self.read_pn_answer}\t{self.read_chrom_answer}\t{self.read_pos_answer}\t{self.read_len_answer}\t{self.read_seq_answer}\t{self.read_quals_answer}\t{self.read_map_q_answer}"
        self.assertEqual(self.dbam_read.to_string(), to_string_answer, f"Both answers should have been {to_string_answer}")
    
    def test_getAsFastQSeq_pairnum(self):
        fastqPnAnswer = f"@{self.read_id_answer}/{self.read_pn_answer}\n{self.read_seq_answer}\n+\n{self.read_quals_answer}\n"
        self.assertEqual(self.dbam_read.get_as_fastq_seq(True), fastqPnAnswer, f"Both answers should have been {fastqPnAnswer}")
    
    def test_getAsFastQSeq_nopairnum(self):
        fastqAnswer = f"@{self.read_id_answer}\n{self.read_seq_answer}\n+\n{self.read_quals_answer}\n"
        self.assertEqual(self.dbam_read.get_as_fastq_seq(), fastqAnswer, f"Both answers should have been {fastqAnswer}")
