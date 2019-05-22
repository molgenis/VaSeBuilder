#!/usr/bin/env python
# Import required libraries/modules
from datetime import datetime
import unittest
from unittest import mock
from unittest.mock import patch, mock_open, call
import os
import sys
import pysam

# Import required class
from VaSeBuilder import VaSeBuilder


# Unittest class for the VaSeBuilder class.
class TestVaSeBuilder(unittest.TestCase):

    # Set up requirements for this test class.
    def setUp(self):
        self.vsBuilder = VaSeBuilder("aap")
        self.varConMap = {'SNP16_247990': ['16', 247986, 56508478]}

    # Tests that the variant type is indeed a SNP
    def test_determineVariantType_snp(self):
        variant_ref = 'C'
        variant_alts = ('G', 'T')
        variant_type_answer = 'snp'
        self.assertEqual(self.vsBuilder.determine_variant_type(variant_ref, variant_alts), variant_type_answer, f"The determined variant type should have been an SNP")

    # Tests that an indel is determined correctly by means of the alternative alleles
    def test_determineVariantType_indel(self):
        variant_ref = 'C'
        variant_alts = ('G', 'TCGATGC')
        variant_type_answer = 'indel'
        self.assertEqual(self.vsBuilder.determine_variant_type(variant_ref, variant_alts), variant_type_answer, f"The determined variant type should have been an indel")

    # Tests that an indel is determined correctly
    def test_determineVariantType_indel2(self):
        variant_ref = 'C,CGATC'
        variant_alts = ('C')
        variant_type_answer = 'indel'
        self.assertEqual(self.vsBuilder.determine_variant_type(variant_ref, variant_alts), variant_type_answer, f"The determined variant type should have been an indel")

    # Tests determining the search window for an SNP
    def test_determineReadSearchWindow(self):
        print("aap")

    # Tests determining the indel read range
    def test_determineIndelReadRange(self):
        variant_position = 200
        variant_ref = 'C'
        variant_alts = ('G', 'TAGCAT')
        indel_range_answer = [variant_position, variant_position+len(variant_alts[1])]
        self.assertListEqual(self.vsBuilder.determine_indel_read_range(variant_position, variant_ref, variant_alts), indel_range_answer, f"The indel read range should have been: {indel_range_answer}")

    # Tests that two variant reads are obtained
    def test_getVariantReads_pos(self):
        read_idlist_answer = ['SRR1039513.12406160', 'SRR1039513.12406160']
        variant_readlist = self.vsBuilder.get_variant_reads("16", 247990, 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
        variant_read_names = [x.get_bam_read_id() for x in variant_readlist]
        self.assertListEqual(variant_read_names, read_idlist_answer, "The lists should be identical but are not")

    # Tests that no reads are obtained
    def test_getVariantReads_neg(self):
        self.assertListEqual(self.vsBuilder.get_variant_reads("16", 1, 1, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb')), [], "Both list should be empty")

    # Tests that the context of BAM reads associated to a variant is determined correctly.
    def test_determineContext_pos(self):
        context_answer = ["16", 247990,  247985, 56508477]
        variant_reads = self.vsBuilder.get_variant_reads("16", 247990, 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
        self.assertListEqual(self.vsBuilder.determine_context(variant_reads, 247990), context_answer, f"The obtained context should have been: {context_answer}.")

    # Test that the context of none existing BAM reads has no context.
    def test_determineContext_neg(self):
        answerList = []
        varReads = self.vsBuilder.get_variant_reads("16", 1, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
        resultList = self.vsBuilder.determine_context(varReads)
        self.assertListEqual(resultList, answerList, "")

    # Tests determining the largest context
    def test_determineLargestContext(self):
        acceptor_context = ['21', 200, 100, 450]
        donor_context = ['21', 200, 150, 500]
        context_answer = ['21', 200, 100, 500]
        self.assertListEqual(self.vsBuilder.determine_largest_context(200, acceptor_context, donor_context), context_answer, f"")

    # Tests that a read is the required read.
    def test_isRequiredRead_pos(self):
        varReads = self.vsBuilder.get_variant_reads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
        if(varReads[1].is_read1):
            self.assertTrue(self.vsBuilder.is_required_read(varReads[1], 'F'), "Should have been true for forward read")
        else:
            self.assertTrue(self.vsBuilder.is_required_read(varReads[0], 'F'), "Should have been true for forward read")

    # Test that a read is not the required read.
    def test_isRequiredRead_neg(self):
        varReads = self.vsBuilder.get_variant_reads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
        if(varReads[0].is_read2):
            self.assertFalse(self.vsBuilder.is_required_read(varReads[0], "F"), "Should have been false for forward read")
        else:
            self.assertFalse(self.vsBuilder.is_required_read(varReads[1], "F"), "Should have been false for forward read")

    # Tests setting the fastq output path for the F (R1) file
    def test_setFastqOutPath(self):
        output_path = '/aap/noot/mies'
        fastq_out_name = ''

    # Tests that the creation id of the VaSeBuilder object is set correctly
    def test_getCreationId(self):
        creationid_answer = 'piet'
        vasebuilder_obj = VaSeBuilder(creationid_answer)
        self.assertEqual(vasebuilder_obj.get_creation_id(), creationid_answer, f"The returned VaSeBuilder identifier should have been: {creationid_answer}")

    # Tests that the creation date of the VaSeBuilder is set correctly.
    def test_getCreationDate(self):
        creation_data_answer = datetime.now().date()
        vasebuilder_obj = VaSeBuilder('aap')
        self.assertEqual(vasebuilder_obj.get_creation_date(), creation_data_answer, f"The returned VaSeBuilder creation data should have been: {creation_data_answer}")

    # Tests that a saved context for a specified variant is indeed returned.
    def test_getVariantContext_pos(self):
        self.vsBuilder.variantContextMap = self.varConMap
        resultList = self.vsBuilder.get_variant_context('SNP16_247990')
        self.vsBuilder.variantContextMap = {}
        self.assertListEqual(resultList, ['16', 247986, 56508478], "Contexts should have been equal")

    # Tests that a None context is returned for a non existent variant.
    def test_getVariantContext_neg(self):
        self.vsBuilder.variantContextMap = self.varConMap
        resultList = self.vsBuilder.get_variant_context('SNP15_10000')
        self.vsBuilder.variantContextMap = {}
        self.assertIsNone(resultList)
