import unittest
from DonorBamRead import DonorBamRead


class TestDonorBamRead(unittest.TestCase):
    # Creates the answer variables and the DonorBamRead used for testing
    def setUp(self):
        self.read_id_answer = "HHKY2CCXX160108:1:2122:24160:2522"
        self.read_pn_answer = "1"
        self.read_chrom_answer = "21"
        self.read_pos_answer = 9411193
        self.read_len_answer = 151
        self.read_seq_answer = "AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGA" \
                               "CTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA"
        self.read_quals_answer = "><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>" \
                                 ">==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>"
        self.read_map_q_answer = 40
        self.dbam_read = DonorBamRead(self.read_id_answer, self.read_pn_answer, self.read_chrom_answer,
                                      self.read_pos_answer, self.read_len_answer, self.read_seq_answer,
                                      self.read_quals_answer, self.read_map_q_answer)

    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_get_bam_read_id(self):
        self.assertEqual(self.dbam_read.get_bam_read_id(), self.read_id_answer, "Both BAM read identifiers should have "
                         f"been {self.read_id_answer}")

    def test_get_bam_read_pair_number(self):
        self.assertEquals(self.dbam_read.get_bam_read_pair_number(), self.read_pn_answer, "Both read pair numbers "
                          f"should have been {self.read_pn_answer}")

    def test_get_bam_read_chrom(self):
        self.assertEquals(self.dbam_read.get_bam_read_chrom(), self.read_chrom_answer, "Both read chromosomes should "
                          f"have been {self.read_chrom_answer}")

    def test_get_bam_read_ref_pos(self):
        self.assertEqual(self.dbam_read.get_bam_read_ref_pos(), self.read_pos_answer, "Both read positions should have "
                         f"been {self.read_pos_answer}")

    def test_get_bam_read_length(self):
        self.assertEqual(self.dbam_read.get_bam_read_length(), self.read_len_answer, "Both read lengths should have "
                         f"been {self.read_len_answer}")

    def test_get_bam_read_ref_end(self):
        read_end_answer = 9411344
        self.assertEqual(self.dbam_read.get_bam_read_ref_end(), read_end_answer, "Both read end positions should have "
                         f"been {read_end_answer}")

    def test_get_bam_read_qual(self):
        self.assertEqual(self.dbam_read.get_bam_read_qual(), self.read_quals_answer, "Both read qualities should have "
                         f"been {self.read_quals_answer}")

    def test_get_bam_read_q_scores(self):
        qscores_answer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28, 27,
                          28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28, 28,
                          28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 27, 28,
                          28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29, 28, 29, 29, 29,
                          29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28, 29, 29, 29, 28, 29,
                          29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29, 30, 30, 30, 29, 30, 30,
                          30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
        self.assertEqual(self.dbam_read.get_bam_read_q_scores(), qscores_answer, "Both Q-scores should have been "
                         f"{qscores_answer}")

    def test_get_bam_read_mapq(self):
        self.assertEqual(self.dbam_read.get_mapping_qual(), self.read_map_q_answer, "Both MapQ values should have been "
                         f"{self.read_map_q_answer}")

    # ====================PERFORM THE TESTS FOR THE STATISTICS METHODS====================
    def test_get_average_q_score(self):
        avg_qscore_answer = 28.490066225165563
        self.assertEqual(self.dbam_read.get_average_qscore(), avg_qscore_answer, "Both average Q-scores should have "
                         f"been {avg_qscore_answer}")

    def test_get_median_q_score(self):
        med_q_score_answer = 28
        self.assertEqual(self.dbam_read.get_median_qscore(), med_q_score_answer, "Both median Q-scores should have "
                         f"been {med_q_score_answer}")

    # ====================PERFORM THE TESTS FOR CHECKING WHETHER THE READ IS R1 OR R2====================
    def test_is_read1(self):
        self.assertTrue(self.dbam_read.is_read1(), "Donor BAM read should have been read 1")

    def test_is_read2(self):
        self.assertFalse(self.dbam_read.is_read2(), "Donor BAM read should not have been read 2")

    # ====================PERFORM THE TESTS FOR RETURNING BAM READ STRING REPRESENTATIONS====================
    def test_to_string(self):
        to_string_answer = f"{self.read_id_answer}\t{self.read_pn_answer}\t{self.read_chrom_answer}\t" \
                           f"{self.read_pos_answer}\t{self.read_len_answer}\t{self.read_seq_answer}\t" \
                           f"{self.read_quals_answer}\t{self.read_map_q_answer}"
        self.assertEqual(self.dbam_read.to_string(), to_string_answer, "Both answers should have been "
                         f"{to_string_answer}")

    def test_get_as_fast_q_seq_pairnum(self):
        fastq_pn_answer = f"@{self.read_id_answer}/{self.read_pn_answer}\n{self.read_seq_answer}\n+\n" \
                          f"{self.read_quals_answer}\n"
        self.assertEqual(self.dbam_read.get_as_fastq_seq(True), fastq_pn_answer, "Both answers should have been "
                         f"{fastq_pn_answer}")

    def test_get_as_fast_q_seq_nopairnum(self):
        fastq_answer = f"@{self.read_id_answer}\n{self.read_seq_answer}\n+\n{self.read_quals_answer}\n"
        self.assertEqual(self.dbam_read.get_as_fastq_seq(), fastq_answer, "Both answers should have been "
                         f"{fastq_answer}")
