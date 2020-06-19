#!/usr/bin/env python
# Import required libraries/modules
from datetime import datetime
import unittest
from unittest import mock
from unittest.mock import patch, mock_open, call
import os
import sys
import random
import pysam

# Import required class
from VaSeBuilder import VaSeBuilder
from VcfVariant import VcfVariant
from DonorBamRead import DonorBamRead


# Unittest class for the VaSeBuilder class.
class TestVaSeBuilder(unittest.TestCase):

    # Set up requirements for this test class.
    def setUp(self):
        self.vs_builder = VaSeBuilder("aap")
        self.var_con_map = {'SNP16_247990': ["16", 247990, 247986, 56508478]}
        self.vcf_variant = VcfVariant("21", 247990, "C", ("G", "T"), ["PASS"], "snp")
        self.bamfile_to_use = "testdata/valbam/SRR1039513.bam"
        self.filter_to_use = ["aap", "noot", "mies"]
        self.vcfvarlist = []

    # ====================PERFORM THE TESTS FOR THE FILTER METHOD====================
    # Tests that a value is indeed in the inclusion filter
    def test_passes_filter_posincl(self):
        val_to_use = "noot"
        self.assertTrue(self.vs_builder.passes_filter(val_to_use, self.filter_to_use), f"Value {val_to_use} should "
                        f"have been in inclusion filter {self.filter_to_use} and therefore return True")

    # Tests that a value is not in the inclusion filter
    def test_passes_filter_negincl(self):
        val_to_use = "piet"
        self.assertFalse(self.vs_builder.passes_filter(val_to_use, self.filter_to_use), f"Value {val_to_use} should "
                         f"not have been in inclusion filter {self.filter_to_use} and therefore return False")

    # Tests that a value is in the exclusion filter
    def test_passes_filter_posexcl(self):
        val_to_use = "mies"
        self.assertFalse(self.vs_builder.passes_filter(val_to_use, self.filter_to_use, is_exclude_filter=True),
                         f"Value {val_to_use} should have been in exclusion filter {self.filter_to_use} "
                         f"and therefore return False")

    # Tests that a value is not in the exclusion filter
    def test_passes_filter_negexcl(self):
        val_to_use = "klaas"
        self.assertTrue(self.vs_builder.passes_filter(val_to_use, self.filter_to_use, is_exclude_filter=True),
                        f"Value {val_to_use} should not have been in exclusion filter {self.filter_to_use} "
                        f"and therefore return True")

    def test_get_sample_vcf_variants_pos(self):
        vcffile = pysam.VariantFile("testdata/vcfDir/SRR1039508.vcf")
        variant_list_answer = [x.to_string() for x in self.vcfvarlist]
        variant_list_answer.sort()

        obtained_variant_list = self.vs_builder.get_sample_vcf_variants(vcffile)
        obtained_list = [x.to_string() for x in obtained_variant_list]
        obtained_list.sort()

        self.assertListEqual(obtained_list, variant_list_answer, "The obtained list of reads should have been "
                             f"{variant_list_answer}")

    # Tests that the variant type is indeed a SNP
    def test_determine_variant_type_snp(self):
        variant_ref = "C"
        variant_alts = ("G", "T")
        variant_type_answer = "snp"
        self.assertEqual(self.vs_builder.determine_variant_type(variant_ref, variant_alts), variant_type_answer,
                         "The determined variant type should have been an SNP")

    # Tests that an indel is determined correctly by means of the alternative alleles
    def test_determine_variant_type_indel(self):
        variant_ref = "C"
        variant_alts = ("G", "TCGATGC")
        variant_type_answer = "indel"
        self.assertEqual(self.vs_builder.determine_variant_type(variant_ref, variant_alts), variant_type_answer,
                         "The determined variant type should have been an indel")

    # Tests that an indel is determined correctly
    def test_determine_variant_type_indel2(self):
        variant_ref = "C,CGATC"
        variant_alts = ("C",)
        variant_type_answer = "indel"
        self.assertEqual(self.vs_builder.determine_variant_type(variant_ref, variant_alts), variant_type_answer,
                         "The determined variant type should have been an indel")

    # Tests determining the search window for an SNP
    def test_determine_read_search_window_snp(self):
        search_window_answer = [247989, 247991]
        obtained_window = self.vs_builder.determine_read_search_window("snp", self.vcf_variant)
        self.assertListEqual(obtained_window, search_window_answer, "Both snp search windows should have been "
                             f"{search_window_answer}")

    # Tests determining the indel read range
    def test_determine_indel_read_range(self):
        variant_position = 200
        variant_ref = "C"
        variant_alts = ("G", "TAGCAT")
        indel_range_answer = [variant_position, variant_position+len(variant_alts[1])]
        self.assertListEqual(self.vs_builder.determine_indel_read_range(variant_position, variant_ref, variant_alts),
                             indel_range_answer, f"The indel read range should have been: {indel_range_answer}")

    # Tests that two variant reads are obtained
    def test_get_variant_reads_pos(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        read_idlist_answer = ["SRR1039513.12406160", "SRR1039513.12406160"]
        variant_readlist = self.vs_builder.get_variant_reads("16_247990", "16", 247989, 247991, bamfile)
        variant_read_names = [x.get_bam_read_id() for x in variant_readlist]
        bamfile.close()
        self.assertListEqual(variant_read_names, read_idlist_answer, "The lists should be identical but are not")

    # Tests that no reads are obtained
    def test_get_variant_reads_neg(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        readlist_answer = []
        obtained_reads = self.vs_builder.get_variant_reads("16_2", "16", 1, 3, bamfile)
        bamfile.close()
        self.assertListEqual(obtained_reads, readlist_answer, "No reads should have been returned")

    #def test_fetch_mate_read(self):
    #def test_filter_variant_reads(self):
    #def test_read_occurence(self):

    # Tests that the context of BAM reads associated to a variant is determined correctly.
    def test_determine_context_pos(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        context_answer = ["16", 247990, 247985, 56508477]
        variant_reads = self.vs_builder.get_variant_reads("16_247990", "16", 247990, 247990, bamfile)
        bamfile.close()
        self.assertListEqual(self.vs_builder.determine_context(variant_reads, 247990, "16"), context_answer,
                             f"The obtained context should have been: {context_answer}.")

    # Test that the context of none existing BAM reads has no context.
    def test_determine_context_neg(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        answer_list = []
        var_reads = self.vs_builder.get_variant_reads("16_2", "16", 1, 100, bamfile)
        result_list = self.vs_builder.determine_context(var_reads, 1, "16")
        bamfile.close()
        self.assertListEqual(result_list, answer_list, "No reads should have been returned")

    # def test_filter_outliers(self):

    # Tests determining the largest context
    def test_determine_largest_context(self):
        acceptor_context = ["21", 200, 100, 450]
        donor_context = ["21", 200, 150, 500]
        context_answer = ["21", 200, 100, 500]
        self.assertListEqual(self.vs_builder.determine_largest_context(200, acceptor_context, donor_context), context_answer, f"")

    # Tests that a read is the required read.
    def test_is_required_read_pos(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        var_reads = self.vs_builder.get_variant_reads("16_247990", "16", 247989, 247991, bamfile)
        bamfile.close()
        if var_reads[1].is_read1:
            self.assertTrue(self.vs_builder.is_required_read(var_reads[1], "F"),
                            "Should have been True for forward read")
        else:
            self.assertTrue(self.vs_builder.is_required_read(var_reads[0], "F"),
                            "Should have been True for forward read")

    # Test that a read is not the required read.
    def test_is_required_read_neg(self):
        bamfile = pysam.AlignmentFile(self.bamfile_to_use, "rb")
        var_reads = self.vs_builder.get_variant_reads("16_247990", "16", 247989, 247991, bamfile)
        bamfile.close()
        if var_reads[0].is_read2:
            self.assertFalse(self.vs_builder.is_required_read(var_reads[0], "F"),
                             "Should have been false for forward read")
        else:
            self.assertFalse(self.vs_builder.is_required_read(var_reads[1], "F"),
                             "Should have been false for forward read")

    # Tests setting the fastq output path for the F (R1) file
    def test_set_fastq_out_path_r1(self):
        outprefix = "/test/outdata/VaSe"
        lanenum = 3
        outpath_answer = f"{outprefix}_{datetime.now().date()}_L{lanenum}_R1.fastq"
        self.assertEqual(self.vs_builder.set_fastq_out_path(outprefix, '1', lanenum), outpath_answer,
                         f"The returned output path for R1 should have been {outpath_answer}")

    def test_set_fastq_out_path_r2(self):
        outprefix = "/test/outdata/VaSe"
        lanenum = 3
        outpath_answer = f"{outprefix}_{datetime.now().date()}_L{lanenum}_R2.fastq"
        self.assertEqual(self.vs_builder.set_fastq_out_path(outprefix, '2', lanenum), outpath_answer,
                         f"The returned output path for R2 should have been {outpath_answer}")

    # Tests that the creation id of the VaSeBuilder object is set correctly
    def test_get_creation_id(self):
        creationid_answer = "piet"
        vasebuilder_obj = VaSeBuilder(creationid_answer)
        self.assertEqual(vasebuilder_obj.get_creation_id(), creationid_answer,
                         f"The returned VaSeBuilder identifier should have been: {creationid_answer}")

    # Tests that a saved context for a specified variant is indeed returned.
    def test_get_variant_context_pos(self):
        self.vs_builder.variantContextMap = self.var_con_map
        result_list = self.vs_builder.get_variant_context("SNP16_247990")
        self.vs_builder.variantContextMap = {}
        self.assertListEqual(result_list, ["16", 247986, 56508478], "Contexts should have been equal")

    # Tests that a None context is returned for a non existent variant.
    def test_get_variant_context_neg(self):
        self.vs_builder.variantContextMap = self.var_con_map
        result_list = self.vs_builder.get_variant_context("SNP15_10000")
        self.vs_builder.variantContextMap = {}
        self.assertIsNone(result_list)

    # =
    # Tests that a list of donor is divided without a remainder list
    def test_divide_donorfastqs_over_acceptors_equal(self):
        donorfq_list = ['aR1.fq', 'bR1.fq', 'cR1.fq', 'dR1.fq', 'eR1.fq', 'fR1.fq']
        divlist_answer = [['fR1.fq', 'cR1.fq'], ['eR1.fq', 'bR1.fq'], ['dR1.fq', 'aR1.fq']]
        self.assertListEqual(self.vs_builder.divide_donorfastqs_over_acceptors(donorfq_list, 3), divlist_answer,
                             f"The returned divided donor fastq list should have been {divlist_answer}")

    # Tests that a list of donors is divided with a remainder list
    def test_divide_donorfastqs_over_acceptors_unequal(self):
        donorfq_list = ['aR1.fq', 'bR1.fq', 'cR1.fq', 'dR1.fq', 'eR1.fq']
        divlist_answer = [['eR1.fq', 'cR1.fq', 'aR1.fq'], ['dR1.fq', 'bR1.fq']]
        self.assertListEqual(self.vs_builder.divide_donorfastqs_over_acceptors(donorfq_list, 2), divlist_answer,
                             f"The returned divided donor fastq list should have been {divlist_answer}")

    # ====================PERFORM TESTS FOR THE SHUFFLING METHODS====================
    # Tests that the shuffling of donor read works as expected
    def test_shuffle_donor_read_identifiers(self):
        read_ids = ["Read4", "Read3", "Read5", "Read1", "Read2"]

        shuffled_answer = ["Read3", "Read2", "Read4", "Read5", "Read1"]
        shuffled_answer = read_ids.copy()
        shuffled_answer.sort()
        random.seed(2)
        random.shuffle(shuffled_answer)

        self.assertListEqual(self.vs_builder.shuffle_donor_read_identifiers(read_ids, 2), shuffled_answer,
                             f"The correct shuffled donor read ids should have been {shuffled_answer}")

    def test_shuffle_donor_add_positions(self):
        add_positions = [6, 8, 2, 0, 1]
        shuffled_answer = add_positions.copy()

    def test_get_saved_insert_position_r1(self):
        read_data = ('1', 100, '2', 100)
        r1pos_answer = 100
        self.assertEqual(self.vs_builder.get_saved_insert_position('1', read_data), r1pos_answer,
                         f"The R1 read insert position should have been {r1pos_answer}")

    def test_get_saved_insert_position_r2(self):
        read_data = ('1', 100, '2', 150)
        r2pos_answer = 150
        self.assertEqual(self.vs_builder.get_saved_insert_position('2', read_data), r2pos_answer,
                         f"The R2 read insert position should have been {r2pos_answer}")

    def test_get_saved_insert_position_nor1(self):
        read_data = ('2', 150)
        r1pos_answer = 'NA'
        self.assertEqual(self.vs_builder.get_saved_insert_position('1', read_data), r1pos_answer,
                         "The R1 read insert position should have been NA")

    def test_get_saved_insert_position_nor2(self):
        read_data = ('1', 100)
        r2pos_answer = 'NA'
        self.assertEqual(self.vs_builder.get_saved_insert_position('2', read_data), r2pos_answer,
                         "The R2 read insert position should have been NA")

    # Adds a donor insert data position for a new fastq
    def test_add_donor_insert_data_addnewfq(self):
        donor_add_data = {}
        fq_name = "ifq.fastq"
        read_id = "iReadId1"
        read_pairnum = "1"
        insert_pos = 100
        donor_add_data_answer = {"ifq.fastq": [("1", 100)]}
        self.vs_builder.add_donor_insert_data(fq_name, read_id, read_pairnum, insert_pos, donor_add_data)
        self.assertDictEqual(donor_add_data, donor_add_data_answer,
                             f"Both donor read add data should have been {donor_add_data_answer}")

    # Tests adding a new donor insert position data for an existing fastq
    def test_add_donor_insert_data_addnewread(self):
        donor_add_data = {"ifq.fastq": []}
        fq_name = "ifq.fastq"
        read_id = "iReadId1"
        read_pairnum = "1"
        insert_pos = 100
        donor_add_data_answer = {"ifq.fastq": [("1", 100)]}
        self.vs_builder.add_donor_insert_data(fq_name, read_id, read_pairnum, insert_pos, donor_add_data)
        self.assertDictEqual(donor_add_data, donor_add_data_answer,
                             f"Both donor read add data should have been {donor_add_data_answer}")

    def test_add_donor_insert_data_addadditional(self):
        donor_add_data = {"ifq.fastq": [("1", 100)]}
        fq_name = "ifq.fastq"
        read_id = "iReadId1"
        read_pairnum = "2"
        insert_pos = 100
        donor_add_data_answer = {"ifq.fastq": [("1", 100, "2", 100)]}
        self.vs_builder.add_donor_insert_data(fq_name, read_id, read_pairnum, insert_pos, donor_add_data)
        self.assertDictEqual(donor_add_data, donor_add_data_answer,
                             f"Both donor read add data should have been {donor_add_data_answer}")

    def test_link_donor_addpos_reads_v2(self):
        addposlist = [6, 8]
        readidlist = ['dRead1', 'dRead2']
        donor_reads = [('dRead1', '1', 'ACAT', '<<>>'), ('dRead1', '2', 'TACA', '>><<'),
                       ('dRead2', '1', 'CAAT', '<><>'), ('dRead2', '2', 'TAAC', '><><')]
        donor_addpos_link_answer = {6: [('dRead1', '1', 'ACAT', '<<>>'), ('dRead1', '2', 'TACA', '>><<')],
                                    8: [('dRead2', '1', 'CAAT', '<><>'), ('dRead2', '2', 'TAAC', '><><')]}
        self.assertDictEqual(self.vs_builder.link_donor_addpos_reads_v2(addposlist, readidlist, donor_reads),
                             donor_addpos_link_answer,
                             f"Both donor read/add position link maps should have been {donor_addpos_link_answer}")

    def test_check_template_size(self):
        num_of_template_reads = 56732
        obtained_template_size = self.vs_builder.check_template_size("testdata/fqDir/SRR1039513_1.fastq.gz")
        self.assertEqual(obtained_template_size, num_of_template_reads,
                         f"The number of template reads was expected to be {num_of_template_reads}")
    
    def test_merge_variants_contexts(self):

    def test_merge_overlap_contexts(self):
        overlap_context_a = None
        overlap_context_b = None
