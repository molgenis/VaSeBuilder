import unittest
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext


class TestOverlapContext(unittest.TestCase):
    # Creates the variables that are needed for each test method
    def setUp(self):
        # Create three DonorBamRead objects to add to the overlap context
        self.read_id_answer = "HHKY2CCXX160108:1:2122:24160:2522"
        self.read_flag_answer = "143"
        self.read_pn_answer = "1"
        self.read_chrom_answer = "21"
        self.read_pos_answer = 9411193
        self.read_len_answer = 151
        self.read_end_answer = 9411331
        self.read_cigarstring_answer = "13S138M"
        self.read_rnext_answer = "21"
        self.read_pnext_answer = 9411300
        self.read_tlen_answer = 4398
        self.read_seq_answer = "AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGA" \
                               "CTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA"
        self.read_quals_answer = "><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>" \
                                 ">==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>"
        self.read_mapq_answer = 40
        self.con_read1 = DonorBamRead(self.read_id_answer, self.read_flag_answer, self.read_pn_answer,
                                      self.read_chrom_answer, self.read_pos_answer, self.read_len_answer,
                                      self.read_end_answer, self.read_cigarstring_answer, self.read_rnext_answer,
                                      self.read_pnext_answer, self.read_tlen_answer, self.read_seq_answer,
                                      self.read_quals_answer, self.read_mapq_answer)
        self.con_read2 = DonorBamRead(self.read_id_answer, self.read_flag_answer, self.read_pn_answer,
                                      self.read_chrom_answer, self.read_pos_answer, self.read_len_answer,
                                      self.read_end_answer, self.read_cigarstring_answer, self.read_rnext_answer,
                                      self.read_pnext_answer, self.read_tlen_answer, self.read_seq_answer,
                                      self.read_quals_answer, self.read_mapq_answer)
        self.con_read3 = DonorBamRead(self.read_id_answer, self.read_flag_answer, self.read_pn_answer,
                                      self.read_chrom_answer, self.read_pos_answer, self.read_len_answer,
                                      self.read_end_answer, self.read_cigarstring_answer, self.read_rnext_answer,
                                      self.read_pnext_answer, self.read_tlen_answer, self.read_seq_answer,
                                      self.read_quals_answer, self.read_mapq_answer)

        # Create the overlap context to test
        self.context_id_answer = '21_9411259'
        self.context_sample_answer = 'testSample'
        self.context_chrom_answer = '21'
        self.context_origin_answer = 9411259
        self.context_start_answer = 9411193
        self.context_end_answer = 9411344
        self.context_len_answer = 151
        self.context_bam_reads_answer = [self.con_read1, self.con_read2, self.con_read3]
        self.unmapped_answer = ["uread_1", "uread_2", "uread_3"]
        self.overlap_context = OverlapContext(self.context_id_answer, self.context_sample_answer,
                                              self.context_chrom_answer, self.context_origin_answer,
                                              self.context_start_answer, self.context_end_answer,
                                              self.context_bam_reads_answer)

        # Create the variables containing the context read test answers
        self.qscore_answer = [29, 27, 28, 30, 30, 30, 29, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 27, 28,
                              27, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28,
                              28, 28, 28, 27, 28, 28, 27, 28, 28, 28, 28, 28, 28, 28, 28, 27, 28, 28, 28, 28, 28, 27,
                              28, 27, 28, 28, 27, 28, 28, 29, 29, 28, 28, 29, 29, 29, 29, 28, 29, 29, 28, 28, 29, 29,
                              28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 29, 29, 29, 28, 29, 29, 29, 28,
                              29, 29, 29, 28, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30, 30, 29, 28, 29, 29,
                              30, 30, 30, 29, 30, 30, 30, 30, 30, 30, 31, 31, 31, 30, 29, 27, 25, 23, 29]
        self.avg_med_read_len_answer = [151, 151]
        self.avg_med_read_qual_answer = [28.490066225165563, 28.490066225165563]
        self.avg_med_read_mapq_answer = [40, 40]

    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_get_context(self):
        context_answer = [self.context_chrom_answer, self.context_origin_answer, self.context_start_answer,
                          self.context_end_answer]
        self.assertListEqual(self.overlap_context.get_context(), context_answer,
                             f"The returned context array should have been {context_answer}")

    def test_get_context_id(self):
        self.assertEqual(self.overlap_context.get_context_id(), self.context_id_answer, "The context id should have "
                         f"been {self.context_id_answer}")

    def test_get_sample_id(self):
        self.assertEqual(self.overlap_context.get_sample_id(), self.context_sample_answer, "The sample id should have "
                         f"been {self.context_sample_answer}")

    def test_get_context_chrom(self):
        self.assertEqual(self.overlap_context.get_context_chrom(), self.context_chrom_answer, "The context chromosome "
                         f"should have been {self.context_chrom_answer}")

    def test_get_context_origin(self):
        self.assertEqual(self.overlap_context.get_context_origin(), self.context_origin_answer, "The context origin "
                         f"should have been {self.context_origin_answer}")

    def test_get_context_start(self):
        self.assertEqual(self.overlap_context.get_context_start(), self.context_start_answer, "The context start "
                         f"position should have been {self.context_start_answer}")

    def test_get_context_end(self):
        self.assertEqual(self.overlap_context.get_context_end(), self.context_end_answer, "Both context end positions "
                         f"should have been {self.context_end_answer}")

    def test_get_context_bam_reads(self):
        self.assertListEqual(self.overlap_context.get_context_bam_reads(), self.context_bam_reads_answer, "The BAM "
                             f"read list should have been {self.context_bam_reads_answer}")

    def test_get_unmapped_read_mate_ids(self):
        self.assertListEqual(self.overlap_context.get_unmapped_read_mate_ids(), self.unmapped_answer, "The unmapped "
                             f"mate read id list should have been {self.unmapped_answer}")

    # ====================PERFORM THE TESTS FOR THE OTHER GETTER METHODS====================
    def test_get_context_length(self):
        self.assertEqual(self.overlap_context.get_context_length(), self.context_len_answer, "The context length "
                         f"should have been {self.context_len_answer}")

    def test_get_start_distance_from_origin(self):
        context_start_distance_answer = self.context_origin_answer - self.context_start_answer
        self.assertEqual(self.overlap_context.get_start_distance_from_origin(), context_start_distance_answer, "The "
                         f"start distance from the context origin should have been {context_start_distance_answer}")

    def test_get_end_distance_from_origin(self):
        context_end_distance_answer = self.context_end_answer - self.context_origin_answer
        self.assertEqual(self.overlap_context.get_end_distance_from_origin(), context_end_distance_answer, "The end "
                         f"distance from the context origin should have been {context_end_distance_answer}")

    # ====================PERFORM THE TESTS FOR GETTING CONTEXT READ INFO====================
    def test_get_number_of_context_reads(self):
        num_of_context_reads_answer = len(self.context_bam_reads_answer)
        self.assertEqual(self.overlap_context.get_number_of_context_reads(), num_of_context_reads_answer, "The number "
                         f"of context reads should have been {num_of_context_reads_answer}")

    def test_get_context_bam_read_ids(self):
        context_read_ids_answer = ["HHKY2CCXX160108:1:2122:24160:2522"]
        self.assertEqual(self.overlap_context.get_context_bam_read_ids(), context_read_ids_answer,
                         f"The list of context read ids should have been {context_read_ids_answer}")

    def test_get_context_bam_read_starts(self):
        context_starts_answer = [9411193, 9411193, 9411193]
        self.assertListEqual(self.overlap_context.get_context_bam_read_starts(), context_starts_answer,
                             f"The list of start positions should have been {context_starts_answer}")

    def test_get_context_bam_read_left_positions(self):
        context_left_pos_answer = [9411193, 9411193, 9411193]
        self.assertEqual(self.overlap_context.get_context_bam_read_left_positions(), context_left_pos_answer,
                         f"The list of left most positions should have been {context_left_pos_answer}")

    def test_get_context_bam_read_ends(self):
        context_ends_answer = [9411344, 9411344, 9411344]
        self.assertListEqual(self.overlap_context.get_context_bam_read_ends(), context_ends_answer,
                             f"The list of end positions should have been {context_ends_answer}")

    def test_get_context_bam_read_right_positions(self):
        context_right_pos_answer = []
        self.assertEqual(self.overlap_context.get_context_bam_read_right_positions(), context_right_pos_answer,
                         f"The list of right most positions should have been {context_right_pos_answer}")

    def test_get_context_bam_read_lengths(self):
        context_read_lens_answer = [151, 151, 151]
        self.assertEqual(self.overlap_context.get_context_bam_read_lengths(), context_read_lens_answer, "The context "
                         f"read lengths should have been {context_read_lens_answer}")

    def test_get_context_bam_read_seqs(self):
        context_read_seqs_answer = [self.read_seq_answer, self.read_seq_answer, self.read_seq_answer]
        self.assertListEqual(self.overlap_context.get_context_bam_read_seqs(), context_read_seqs_answer, "The context "
                             f"read sequences should have been {context_read_seqs_answer}")

    def test_get_context_bam_read_qualities(self):
        context_read_quals_answer = [self.read_quals_answer, self.read_quals_answer, self.read_quals_answer]
        self.assertEqual(self.overlap_context.get_context_bam_read_qualities(), context_read_quals_answer, "The read "
                         f"qualites should have been {context_read_quals_answer}")

    def test_get_context_bam_read_q_scores(self):
        context_read_q_scores_answer = [self.qscore_answer] + [self.qscore_answer] + [self.qscore_answer]
        self.assertListEqual(self.overlap_context.get_context_bam_read_q_scores(), context_read_q_scores_answer,
                             f"The Q-scores should have been {context_read_q_scores_answer}")

    def test_get_context_bam_read_map_qs(self):
        context_read_map_qs_answer = [40, 40, 40]
        self.assertListEqual(self.overlap_context.get_context_bam_read_map_qs(), context_read_map_qs_answer, "The "
                             f"qualities should have been {context_read_map_qs_answer}")

    def test_read_is_in_context_pos(self):
        self.assertTrue(self.overlap_context.read_is_in_context(self.read_id_answer), f"Read id {self.read_id_answer} "
                        "should have been in the context")
    
    def test_read_is_in_context_neg(self):
        read_in_context_answer_neg = "HHKY2CCXX160108:1:2122:24160:2555"
        self.assertFalse(self.overlap_context.read_is_in_context(read_in_context_answer_neg), "Read id "
                         f"{read_in_context_answer_neg} should not have been in the context")

    # ====================PERFORM THE TESTS FOR UNMAPPED MATE IDS====================
    def test_add_unmapped_mate_id(self):
        read_id_to_add = "uReadId1"
        self.overlap_context.add_unmapped_mate_id(read_id_to_add)
        self.assertTrue(read_id_to_add in self.overlap_context.unmapped_read_mate_ids,
                        f"Read id {read_id_to_add} should have been in the list of reads with unmapped mates")

    def test_set_unmapped_mate_ids(self):
        unmapped_read_ids = ["uReadId1", "uReadId2", "uReadId3"]
        self.overlap_context.set_unmapped_mate_ids(unmapped_read_ids)
        self.assertListEqual(self.overlap_context.unmapped_read_mate_ids, unmapped_read_ids,
                             f"The list of unmapped read mate ids should have been {unmapped_read_ids}")

    def test_read_has_unmapped_mate_true(self):
        read_id_to_search = "uread_1"
        self.assertTrue(self.overlap_context.read_has_unmapped_mate(read_id_to_search),
                        f"Read {read_id_to_search} should have an unmapped mate")

    def test_read_has_unmapped_mate_false(self):
        read_id_to_search = "fReadId1"
        self.assertFalse(self.overlap_context.read_has_unmapped_mate(read_id_to_search),
                         f"Read {read_id_to_search} should not have an unmapped mate")

    # ====================PERFORM THE TESTS FOR THE OVERLAP CONTEXT STATISTICS====================
    def test_get_average_and_median_read_length(self):
        self.assertListEqual(self.overlap_context.get_average_and_median_read_length(), self.avg_med_read_len_answer,
                             f"Average and median read length should have been {self.avg_med_read_len_answer}")

    def test_get_average_and_median_read_qual(self):
        self.assertListEqual(self.overlap_context.get_average_and_median_read_qual(), self.avg_med_read_qual_answer,
                             f"Average and median read length should have been {self.avg_med_read_qual_answer}")

    def test_get_average_and_median_read_map_q(self):
        self.assertListEqual(self.overlap_context.get_average_and_median_read_map_q(), self.avg_med_read_mapq_answer,
                             f"Average and median read length should have been {self.avg_med_read_mapq_answer}")

    # ====================PERFOM THE TESTS FOR THE OTHER METHODS====================
    def test_to_string(self):
        to_string_answer = f"{self.context_id_answer}\t{self.context_sample_answer}\t{self.context_chrom_answer}\t" \
                           f"{self.context_origin_answer}\t{self.context_start_answer}\t{self.context_end_answer}\t" \
                           f"{len(self.context_bam_reads_answer)}\t" + \
                           ";".join([x.get_bam_read_id() for x in self.context_bam_reads_answer])
        self.assertEqual(self.overlap_context.to_string(), to_string_answer, "String representation should have "
                         f"been {to_string_answer}")

    def test_to_statistics_string(self):
        to_statistics_answer = f"{self.context_id_answer}\t{self.avg_med_read_len_answer[0]}\t" \
                               f"{self.avg_med_read_len_answer[1]}\t{self.avg_med_read_qual_answer[0]}\t" \
                               f"{self.avg_med_read_qual_answer[1]}\t{self.avg_med_read_mapq_answer[0]}\t" \
                               f"{self.avg_med_read_mapq_answer[1]}"
        self.assertEqual(self.overlap_context.to_statistics_string(), to_statistics_answer, "Statistics string should "
                         f"have been {to_statistics_answer}")

    def test_compare(self):
        compare_answer = {}
        self.assertDictEqual(self.overlap_context.compare(self.overlap_context), compare_answer, "Compare results"
                             f"dictionary should have been {compare_answer}")
