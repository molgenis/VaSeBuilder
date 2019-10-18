import unittest
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext


class TestVariantContext(unittest.TestCase):
    def setUp(self):
        # Create the variables saving the context BAM reads
        self.read_id_answer = "HHKY2CCXX160108:1:2122:24160:2522"
        self.read_donor_id_answer = "HHKY2CCXX160108:1:2122:24160:2555"
        self.read_flag_answer = "143"
        self.read_pn_answer = "1"
        self.read_donor_pn_answer = "2"
        self.read_chrom_answer = "21"
        self.read_pos_answer = 9411193
        self.read_len_answer = 151
        self.read_donor_len_answer = 108
        self.read_end_answer = 9411334
        self.read_cigarstring_answer = "13S138M"
        self.read_donor_cigarstring_answer = "43H108M"
        self.read_rnext_answer = "21"
        self.read_pnext_answer = 9511300
        self.read_tlen_answer = 400
        self.read_seq_answer = "AGAAAAAGTCTTTAGATGGGATCTTCCTCCAAAGAAATTGTAGTTTTCTTCTGGCTTAGAGGTAGATCATCTTGGTCCAATCAGA" \
                               "CTGAAATGCCTTGAGGCTAGATTTCAGTCTTTGTGGCAGCTGGTGAATTTCTAGTTTGCCTTTTCA"
        self.read_quals_answer = "><=???>==<=====<====<=<==<==<=====<============<==<========<=====<=<==<==>>==>>>>=>" \
                                 ">==>>=>>>>>>>>>=>>>>>>=>>>=>>>=>>>>>?????????>=>>???>??????@@@?><:8>"
        self.read_map_q_answer = 40
        self.acceptor_read = DonorBamRead(self.read_id_answer, self.read_flag_answer, self.read_pn_answer,
                                          self.read_chrom_answer, self.read_pos_answer, self.read_len_answer,
                                          self.read_end_answer, self.read_cigarstring_answer, self.read_rnext_answer,
                                          self.read_pnext_answer, self.read_tlen_answer, self.read_seq_answer,
                                          self.read_quals_answer, self.read_map_q_answer)
        self.donor_read = DonorBamRead(self.read_donor_id_answer, self.read_flag_answer, self.read_donor_pn_answer,
                                       self.read_chrom_answer, self.read_pos_answer, self.read_donor_len_answer,
                                       self.read_end_answer, self.read_donor_cigarstring_answer, self.read_rnext_answer,
                                       self.read_pnext_answer, self.read_tlen_answer, self.read_seq_answer,
                                       self.read_quals_answer, self.read_map_q_answer)

        self.acceptor_read_ids_answer = ["HHKY2CCXX160108:1:2122:24160:2522"]
        self.donor_read_ids_answer = ["HHKY2CCXX160108:1:2122:24160:2555"]
        self.read_lens_answer = [self.read_len_answer, self.read_len_answer, self.read_len_answer]
        self.donor_read_lens_answer = [self.read_donor_len_answer, self.read_donor_len_answer,
                                       self.read_donor_len_answer]

        # Create the variables containing the context answers
        self.context_id_answer = "21_9411259"
        self.context_sample_answer = "testsample"
        self.context_chrom_answer = "21"
        self.context_origin_answer = 9411259
        self.acceptor_start_answer = 9411210
        self.donor_start_answer = 9411193
        self.context_start_answer = 9411193
        self.acceptor_end_answer = 9411344
        self.donor_end_answer = 9411301
        self.context_end_answer = 9411344
        self.avg_med_len_answer = []
        self.avg_med_qual_answer = []
        self.avg_med_map_q_answer = []

        # Creat the variables containing the acceptor context answers
        self.acceptor_reads_answer = [self.acceptor_read, self.acceptor_read, self.acceptor_read]
        self.acceptor_read_starts_answer = [self.read_pos_answer, self.read_pos_answer, self.read_pos_answer]
        self.acceptor_read_ends_answer = [9411344, 9411344, 9411344]
        self.acceptor_right_pos_answer = []

        # Create the variables containing the donor context answers
        self.donor_reads_answer = [self.donor_read, self.donor_read, self.donor_read]
        self.donor_read_starts_answer = [9411193, 9411193, 9411193]
        self.donor_read_ends_answer = [9411301, 9411301, 9411301]
        self.donor_left_pos_answer = []

        # Create the variables containing the accept, donor and variant context
        self.acceptor_context_answer = OverlapContext(self.context_id_answer, self.context_sample_answer,
                                                      self.context_chrom_answer, self.context_origin_answer,
                                                      self.acceptor_start_answer, self.acceptor_end_answer,
                                                      self.acceptor_reads_answer)
        self.donor_context_answer = OverlapContext(self.context_id_answer, self.context_sample_answer,
                                                   self.context_chrom_answer, self.context_origin_answer,
                                                   self.donor_start_answer, self.donor_end_answer,
                                                   self.donor_reads_answer)
        self.variant_context = VariantContext(self.context_id_answer, self.context_sample_answer,
                                              self.context_chrom_answer, self.context_origin_answer,
                                              self.context_start_answer, self.context_end_answer,
                                              self.acceptor_reads_answer, self.donor_reads_answer,
                                              self.acceptor_context_answer, self.donor_context_answer)

        # Create the variables containing the string output answers
        self.to_statistics_answer = f"{self.context_id_answer}\t"

    # ====================PERFORM THE TESTS FOR THE GETTER METHODS====================
    def test_get_variant_context_id(self):
        self.assertEqual(self.variant_context.get_variant_context_id(), self.context_id_answer, "The variant context id"
                         f" should have been {self.context_id_answer}")

    def test_get_variant_context_sample(self):
        self.assertEqual(self.variant_context.get_variant_context_sample(), self.context_sample_answer, "The variant "
                         f"context sample should have been {self.context_sample_answer}")

    def test_get_variant_context_chrom(self):
        self.assertEqual(self.variant_context.get_variant_context_chrom(), self.context_chrom_answer, "The variant "
                         f"context chromosome should have been {self.context_chrom_answer}")

    def test_get_variant_context_origin(self):
        self.assertEqual(self.variant_context.get_variant_context_origin(), self.context_origin_answer, "The variant "
                         f"context origin position should have been {self.context_origin_answer}")

    def test_get_variant_context_start(self):
        self.assertEqual(self.variant_context.get_variant_context_start(), self.context_start_answer, "The variant "
                         f"context start position should have been {self.context_start_answer}")

    def test_get_variant_context_end(self):
        self.assertEqual(self.variant_context.get_variant_context_end(), self.context_end_answer, "The variant context "
                         f"end position should have been {self.context_end_answer}")

    def test_get_variant_context_acceptor_reads(self):
        self.assertListEqual(self.variant_context.get_acceptor_reads(), self.acceptor_reads_answer, "The variant "
                                                                                                    "context acceptor "
                                                                                                    "reads are not what"
                                                                                                    " was expected")
    
    def test_get_variant_context_donor_reads(self):
        self.assertListEqual(self.variant_context.get_donor_reads(), self.donor_reads_answer, "The variant context "
                                                                                              "donor reads are not what"
                                                                                              " was expected")

    def test_get_acceptor_context(self):
        self.assertEqual(self.variant_context.get_acceptor_context(), self.acceptor_context_answer, "The acceptor "
                                                                                                    "context is not "
                                                                                                    "what was expected")

    def test_get_donor_context(self):
        self.assertEqual(self.variant_context.get_donor_context(), self.donor_context_answer, "The donor context is not"
                                                                                              " what was expected")

    def test_get_unmapped_acceptor_mate_ids(self):
        vua_read_ids = ['vuaRead1', 'vuaRead2', 'vuaRead3']
        self.variant_context.set_unmapped_acceptor_mate_ids(vua_read_ids)
        self.assertListEqual(self.variant_context.get_unmapped_acceptor_mate_ids(), vua_read_ids, f"The variant context"
                             f" unmapped acceptor mate ids should have been {vua_read_ids}")

    def test_get_unmapped_donor_mate_ids(self):
        vud_read_ids = ['vudRead1', 'vudRead2', 'vudRead3']
        self.variant_context.set_unmapped_donor_mate_ids(vud_read_ids)
        self.assertListEqual(self.variant_context.get_unmapped_donor_mate_ids(), vud_read_ids, "The variant context "
                             f"unmapped donor mate ids should have been {vud_read_ids}")

    # ====================PERFORM THE TESTS FOR THE OTHER GETTER METHODS====================
    def test_get_variant_context_length(self):
        self.assertEqual(self.variant_context.get_variant_context_length(),
                         (abs(self.context_end_answer - self.context_start_answer)),
                         "The length of the variant context should have been "
                         f"{abs(self.context_end_answer - self.context_start_answer)}")
    
    def test_get_start_distance_from_origin(self):
        self.assertEqual(self.variant_context.get_start_distance_from_origin(),
                         (abs(self.context_origin_answer - self.context_start_answer)),
                         "The start distance from the origin should have been "
                         f"{abs(self.context_origin_answer - self.context_start_answer)}")
    
    def test_get_end_distance_from_origin(self):
        self.assertEqual(self.variant_context.get_end_distance_from_origin(),
                         (abs(self.context_end_answer - self.context_origin_answer)),
                         "The end distance from the origin should have been "
                         f"{abs(self.context_end_answer - self.context_origin_answer)}")

    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT ACCEPTOR READ DATA====================
    def test_get_number_of_variant_context_acceptor_reads(self):
        self.assertEqual(self.variant_context.get_number_of_acceptor_reads(),
                         len(self.acceptor_reads_answer), "The variant context acceptor reads are not what was "
                                                          "expected")
    
    def test_get_variant_context_acceptor_read_ids(self):
        self.assertListEqual(self.variant_context.get_acceptor_read_ids(), self.acceptor_read_ids_answer, "The variant "
                             f"context acceptor read ids should have been {self.acceptor_read_ids_answer}")
    
    def test_get_variant_context_acceptor_read_starts(self):
        self.assertListEqual(self.variant_context.get_acceptor_read_starts(), self.acceptor_read_starts_answer, "The "
                             f"variant context acceptor read starts should have been "
                             f"{self.acceptor_read_starts_answer}")

    def test_get_acceptor_read_left_positions(self):
        self.assertListEqual(self.variant_context.get_acceptor_read_left_positions(), self.acceptor_read_starts_answer,
                             f"The variant context acceptor read left most positions should have been "
                             f"{self.acceptor_read_starts_answer}")

    def test_get_acceptor_read_ends(self):
        self.assertListEqual(self.variant_context.get_acceptor_read_ends(), self.acceptor_read_ends_answer,
                             f"The variant context acceptor read end positions should have been "
                             f"{self.acceptor_read_ends_answer}")

    def test_get_acceptor_read_right_positions(self):
        self.assertListEqual(self.variant_context.get_acceptor_read_right_positions(), self.acceptor_right_pos_answer,
                             f"The variant context acceptor read right most positions should have been "
                             f"{self.acceptor_right_pos_answer}")

    # ====================PERFORM THE TESTS FOR GETTING THE VARIANT CONTEXT DONOR READ DATA====================
    def test_get_number_of_variant_context_donor_reads(self):
        self.assertEqual(self.variant_context.get_number_of_donor_reads(), len(self.donor_reads_answer),
                         "The variant context donor reads are not what was expected")
    
    def test_get_variant_context_donor_read_ids(self):
        self.assertListEqual(self.variant_context.get_donor_read_ids(), self.donor_read_ids_answer,
                             f"the variant context donor read ids should have been {self.donor_read_ids_answer}")
    
    def test_get_variant_context_donor_read_starts(self):
        self.assertListEqual(self.variant_context.get_donor_read_starts(), self.donor_read_starts_answer,
                             "The variant context donor read start positions should have been "
                             f"{self.donor_read_starts_answer}")
    
    def test_get_variant_context_donor_read_left_positions(self):
        self.assertListEqual(self.variant_context.get_donor_read_left_positions(), self.donor_left_pos_answer,
                             "The variant context donor read left most positions should have been "
                             f"{self.donor_left_pos_answer}")
    
    def test_get_variant_context_donor_read_ends(self):
        self.assertListEqual(self.variant_context.get_donor_read_ends(), self.donor_read_ends_answer,
                             "The variant context donor read end positions should have been "
                             f"{self.donor_read_ends_answer}")
    
    def test_get_variant_context_donor_read_right_positions(self):
        self.assertListEqual(self.variant_context.get_donor_read_right_positions(), self.donor_read_ends_answer,
                             "The variant context donor read right most positions should have been "
                             f"{self.donor_read_ends_answer}")

    # ====================PERFORM THE TESTS FOR ADDING DATA TO THE VARIANT CONTEXT====================
    def test_set_acceptor_context(self):
        acc_con_to_add = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptor_reads_answer)
        self.variant_context.set_acceptor_context(acc_con_to_add)
        self.assertEqual(self.variant_context.get_acceptor_context().to_string(), acc_con_to_add.to_string(),
                         "The saved acceptor context is not the same as the one which was just set")

    def test_set_donor_context(self):
        don_con_to_add = OverlapContext("22_10001", "settest", "22", 10001, 9990, 10100, self.donor_reads_answer)
        self.variant_context.set_donor_context(don_con_to_add)
        self.assertEqual(self.variant_context.get_donor_context().to_string(), don_con_to_add.to_string(),
                         "The saved donor context is not the same as the one which was just set")

    def test_add_acceptor_context(self):
        acc_con_to_add = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.acceptor_reads_answer)
        self.variant_context.add_acceptor_context("22_10001", "settest", "22", 10001, 9990, 10100,
                                                  self.acceptor_reads_answer)
        self.assertEqual(self.variant_context.get_acceptor_context().to_string(), acc_con_to_add.to_string(),
                         "The saved donor context is not the same as the one which was just set")

    def test_add_donor_context(self):
        don_con_to_add = OverlapContext('22_10001', 'settest', '22', 10001, 9990, 10100, self.donor_reads_answer)
        self.variant_context.add_donor_context('22_10001', 'settest', '22', 10001, 9990, 10100, self.donor_reads_answer)
        self.assertEqual(self.variant_context.get_donor_context().to_string(), don_con_to_add.to_string(),
                         "The saved donor context is not the same as the one which was just set")

    # ====================PERFORM THE TESTS FOR GETTING ACCEPTOR CONTEXT DATA====================
    def test_get_acceptor_context_id(self):
        self.assertEqual(self.variant_context.get_acceptor_context_id(), self.context_id_answer, "The acceptor context "
                                                                                                 "id should have been")

    def test_get_acceptor_sample_id(self):
        self.assertEqual(self.variant_context.get_acceptor_context_sample_id(), self.context_sample_answer,
                         f"The acceptor context sample id should have been {self.context_sample_answer}")

    def test_get_acceptor_context_chrom(self):
        self.assertEqual(self.variant_context.get_acceptor_context_chrom(), self.context_chrom_answer,
                         f"The acceptor context chromosome should have been {self.context_chrom_answer}")

    def test_get_acceptor_context_origin(self):
        self.assertEqual(self.variant_context.get_acceptor_context_origin(), self.context_origin_answer,
                         f"The acceptor context origin should have been {self.context_origin_answer}")

    def test_get_acceptor_context_start(self):
        self.assertEqual(self.variant_context.get_acceptor_context_start(), self.acceptor_start_answer,
                         f"The acceptor context start positions should have been {self.acceptor_start_answer}")

    def test_get_acceptor_context_end(self):
        self.assertEqual(self.variant_context.get_acceptor_context_end(), self.acceptor_end_answer,
                         f"The acceptor context end positions should have been {self.acceptor_end_answer}")

    def test_get_acceptor_context_length(self):
        self.assertEqual(self.variant_context.get_acceptor_context_length(), abs(self.acceptor_end_answer -
                                                                                 self.acceptor_start_answer),
                         "The acceptor context length should have been "
                         f"{abs(self.acceptor_end_answer - self.acceptor_start_answer)}")

    def test_get_acceptor_context_reads(self):
        self.assertEqual(self.variant_context.get_acceptor_context_reads(), self.acceptor_reads_answer,
                         "The acceptor context reads are not what was expected")

    def test_get_acceptor_context_read_ids(self):
        self.assertListEqual(self.variant_context.get_acceptor_context_read_ids(), self.acceptor_read_ids_answer,
                             f"The acceptor context read ids should have been {self.acceptor_read_ids_answer}")

    def test_get_acceptor_context_read_starts(self):
        self.assertEqual(self.variant_context.get_acceptor_context_read_starts(), self.acceptor_read_starts_answer,
                         "The acceptor context read start positions should have been "
                         f"{self.acceptor_read_starts_answer}")

    def test_get_acceptor_context_read_left_positions(self):
        self.assertEqual(self.variant_context.get_acceptor_context_read_left_positions(),
                         self.acceptor_read_starts_answer,
                         "The acceptor context read left most positions should have been "
                         f"{self.acceptor_read_starts_answer}")

    def test_get_acceptor_context_read_ends(self):
        self.assertEqual(self.variant_context.get_acceptor_context_read_ends(), self.acceptor_read_ends_answer,
                         f"The acceptor read end positions should have been {self.acceptor_read_ends_answer}")

    def test_get_acceptor_context_read_right_positions(self):
        self.assertEqual(self.variant_context.get_acceptor_context_read_right_positions(),
                         self.acceptor_right_pos_answer, "The acceptor context read right most positions should have "
                         f"been {self.acceptor_right_pos_answer}")

    def test_get_acceptor_context_read_lengths(self):
        self.assertEqual(self.variant_context.get_acceptor_context_read_lengths(), self.read_lens_answer,
                         f"The acceptor context read lengths should have been {self.read_lens_answer}")

    def test_get_acceptor_context_unmapped_mate_ids(self):
        au_read_ids = ['auRead1', 'auRead2', 'auRead3']
        self.variant_context.set_acceptor_context_unmapped_mates(au_read_ids)
        self.assertListEqual(self.variant_context.get_acceptor_context_unmapped_mate_ids(), au_read_ids,
                             f"The acceptor context unmapped ids should have been {au_read_ids}")

    # ====================PERFORM THE TESTS FOR GETTING DONOR CONTEXT DATA====================
    def test_get_donor_context_id(self):
        self.assertEqual(self.variant_context.get_donor_context_id(), self.context_id_answer,
                         f"The donor context id should have been {self.context_id_answer}")

    def test_get_donor_sample_id(self):
        self.assertEqual(self.variant_context.get_donor_context_sample_id(), self.context_sample_answer,
                         f"The donor context sample id should have been {self.context_sample_answer}")

    def test_get_donor_context_chrom(self):
        self.assertEqual(self.variant_context.get_donor_context_chrom(), self.context_chrom_answer,
                         f"The donor context chromosome should have been {self.context_chrom_answer}")

    def test_get_donor_context_origin(self):
        self.assertEqual(self.variant_context.get_donor_context_origin(), self.context_origin_answer,
                         f"The donor context origin should have been {self.context_origin_answer}")

    def test_get_donor_context_start(self):
        self.assertEqual(self.variant_context.get_donor_context_start(), self.donor_start_answer,
                         f"The donor context start position should have been {self.donor_start_answer}")

    def test_get_donor_context_end(self):
        self.assertEqual(self.variant_context.get_donor_context_end(), self.donor_end_answer,
                         f"The donor context end position should have been {self.donor_end_answer}")

    def test_get_donor_context_length(self):
        self.assertEqual(self.variant_context.get_donor_context_length(),
                         abs(self.donor_end_answer - self.donor_start_answer),
                         "The donor context length should have been "
                         f"{abs(self.donor_end_answer - self.donor_start_answer)}")

    def test_get_donor_context_reads(self):
        self.assertEqual(self.variant_context.get_donor_context_reads(), self.donor_reads_answer,
                         "The donor context reads are not what was expected")

    def test_get_donor_context_read_ids(self):
        self.assertListEqual(self.variant_context.get_donor_context_read_ids(), self.donor_read_ids_answer,
                             f"The donor context read ids should have been {self.donor_read_ids_answer}")

    def test_get_donor_context_read_starts(self):
        self.assertEqual(self.variant_context.get_donor_context_read_starts(), self.donor_read_starts_answer,
                         f"The donor context read start positions should have been {self.donor_read_starts_answer}")

    def test_get_donor_context_read_left_positions(self):
        self.assertEqual(self.variant_context.get_donor_context_read_left_positions(), self.donor_left_pos_answer,
                         f"The donor context read left most positions should have been {self.donor_left_pos_answer}")

    def test_get_donor_context_read_ends(self):
        self.assertEqual(self.variant_context.get_donor_context_read_ends(), self.donor_read_ends_answer,
                         f"The donor context read end positions should have been {self.donor_read_ends_answer}")

    def test_get_donor_context_read_right_positions(self):
        self.assertEqual(self.variant_context.get_donor_context_read_right_positions(), self.donor_read_ends_answer,
                         f"The donor context read right most positions should have been {self.donor_read_ends_answer}")

    def test_get_donor_context_read_lengths(self):
        self.assertEqual(self.variant_context.get_donor_context_read_lengths(), self.donor_read_lens_answer,
                         f"The donor context read lengths should have been {self.read_lens_answer}")

    def test_get_donor_context_unmapped_mate_ids(self):
        du_read_ids = ['duRead1', 'duRead2', 'duRead3']
        self.variant_context.set_donor_context_unmapped_mates(du_read_ids)
        self.assertListEqual(self.variant_context.get_donor_context_unmapped_mate_ids(), du_read_ids,
                             f"The donor context id should have been {du_read_ids}")

    # ====================PERFORM THE TESTS FOR THE OUTPUT METHODS====================
    def test_to_string(self):
        to_string_answer = f"{self.context_id_answer}\t{self.context_sample_answer}\t{self.context_chrom_answer}\t" \
                           f"{self.context_origin_answer}\t{self.context_start_answer}\t{self.context_end_answer}\t" \
                           f"{abs(self.acceptor_end_answer - self.acceptor_start_answer)}\t" \
                           f"{abs(self.donor_end_answer - self.donor_start_answer)}\t" \
                           f"{len(self.acceptor_reads_answer)}\t{len(self.donor_reads_answer)}\t" \
                           f"{float(len(self.acceptor_reads_answer) / len(self.donor_reads_answer))}\t" + \
                           ";".join(self.acceptor_read_ids_answer) + "\t" + ';'.join(self.donor_read_ids_answer)
        self.assertEqual(self.variant_context.to_string(), to_string_answer,
                         f"The toString line should have been {to_string_answer}")

    def test_to_statistics_string(self):
        to_statistics_answer = f"{self.context_id_answer}\t"
        self.assertEqual(self.variant_context.to_statistics_string(), self.to_statistics_answer,
                         f"The statistics line should have been {self.to_statistics_answer}")
