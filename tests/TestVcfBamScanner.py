#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os
import sys
import pysam

# Import required class
# sys.path.insert(0,'./../')
from VcfBamScanner import VcfBamScanner


# Unittest class for the VcfBamScanner class
class TestVcfBamScanner(unittest.TestCase):
    
    # Set up some things for the testing methods.
    def setUp(self):
        self.vb_scanner = VcfBamScanner()
        self.vcf_folders = ['testdata/vcfDir']
        self.bam_folders = ['testdata/bamDir']
        self.none_folders = ['testdata/doesnotexist']
        self.valid_bam_file = "testdata/bamDir/SRR1039508.bam"
        self.valid_vcf_file = "testdata/vcfDir/SRR1039508.vcf"

    # Tests the scanning of VCF folders that exist and contain vcf/vcf.gz files
    def test_scan_vcf_folders_pos(self):
        answer_dict = {'SRR1039508': 'testdata/vcfDir/SRR1039508.vcf', 'SRR1039512': 'testdata/vcfDir/SRR1039512.vcf'}
        result_dict = self.vb_scanner.scan_vcf_folders(self.vcf_folders)
        self.assertDictEqual(result_dict, answer_dict, "Dicts should have been the same")

    # Tests the scanning of VCF folders that do not exist / contain no vcf/vcf.gz files.
    def test_scan_vcf_folders_neg(self):
        answer_dict = {}
        result_dict = self.vb_scanner.scan_vcf_folders(self.none_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been empty but are not")

    # Tests the scanning of folders containg BAM files.
    def test_scan_bam_folders_pos(self):
        answer_dict = {'SRR1039508': 'testdata/bamDir/SRR1039508.bam', 'SRR1039512': 'testdata/bamDir/SRR1039512.bam'}
        result_dict = self.vb_scanner.scan_bam_folders(self.bam_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been the same")

    # Tests the scanning of folders that do not exist / contain no BAM files.
    def test_scan_bam_folders_neg(self):
        answer_dict = {}
        result_dict = self.vb_scanner.scan_bam_folders(self.none_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been empty")

    # Tests that a VCF file has sample information
    def test_vcf_has_sample_name_pos(self):
        valid_vcf_file = pysam.VariantFile("testdata/vcfDir/SRR1039508.vcf")
        result_bool = self.vb_scanner.vcf_has_sample_name(valid_vcf_file)
        valid_vcf_file.close()
        self.assertTrue(result_bool, f"VCF file {valid_vcf_file} should have had a sample identifier")

    # Tests that a BAM file has sample information.
    def test_bam_has_sample_name_pos(self):
        bam_file = pysam.AlignmentFile("testdata/bamDir/SRR1039508.bam")
        result_bool = self.vb_scanner.bam_has_sample_name(bam_file)
        bam_file.close()
        self.assertTrue(result_bool)

    # Tests that a BAM file has no sample information.
    def test_bam_has_sample_name_neg(self):
        bam_file = pysam.AlignmentFile("testdata/noSampleDir/noSampleBam.bam")
        result_bool = self.vb_scanner.bam_has_sample_name(bam_file)
        bam_file.close()
        self.assertFalse(result_bool)

    # Tests tha obtaining the BAM sample name is obtained
    def test_get_bam_sample_name_pos(self):
        valid_bam_file = "testdata/bamDir/SRR1039508.bam"
        sample_id_answer = "SRR1039508"
        self.assertEqual(self.vb_scanner.get_bam_sample_name(valid_bam_file), sample_id_answer,
                         f"Both BAM sample ids should have been {sample_id_answer}")

    # Tests that no BAM sample id is obtained
    def test_get_bam_sample_name_neg(self):
        no_sample_bam = "testdata/noSampleDir/noSampleBam.bam"
        self.assertIsNone(self.vb_scanner.get_bam_sample_name(no_sample_bam),
                          f"BAM file {no_sample_bam} should have no sample")

    # Tests whether the VCF to BAM map will be constructed properly.
    def test_get_vcf_to_bam_map(self):
        answer_dict = {'testdata/vcfDir/SRR1039508.vcf': 'testdata/bamDir/SRR1039508.bam',
                       'testdata/vcfDir/SRR1039512.vcf': 'testdata/bamDir/SRR1039512.bam'}
        vcf_fs = self.vb_scanner.scan_vcf_folders(self.vcf_folders)    # Provide VcfBamScanner a valid list of VCF files
        bam_fs = self.vb_scanner.scan_bam_folders(self.bam_folders)    # Provide VcfBamScanner a valid list of BAM files
        result_dict = self.vb_scanner.get_vcf_to_bam_map()
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been the same")

    # Tests that the correct sequence names are obtained from the testfile
    def test_get_alignment_sequence_names(self):
        sequence_names_answer = {"16"}
        bamfile = pysam.AlignmentFile(self.valid_bam_file)
        obtained_sequence_names = self.vb_scanner.get_alignment_sequence_names(bamfile)
        bamfile.close()
        self.assertSetEqual(obtained_sequence_names, sequence_names_answer, "The set of sequence names should have been"
                            f" {sequence_names_answer}")

    # Tests that the contig names from a vcf are obtained correctly
    def test_get_variant_sequence_names(self):
        sequence_names_answer = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16",
                                 "17", "18", "19", "20", "21", "22", "X", "Y", "MT", "GL000207.1", "GL000226.1",
                                 "GL000229.1", "GL000231.1", "GL000210.1", "GL000239.1", "GL000235.1", "GL000201.1",
                                 "GL000247.1", "GL000245.1", "GL000197.1", "GL000203.1", "GL000246.1", "GL000249.1",
                                 "GL000196.1", "GL000248.1", "GL000244.1", "GL000238.1", "GL000202.1", "GL000234.1",
                                 "GL000232.1", "GL000206.1", "GL000240.1", "GL000236.1", "GL000241.1", "GL000243.1",
                                 "GL000242.1", "GL000230.1", "GL000237.1", "GL000233.1", "GL000204.1", "GL000198.1",
                                 "GL000208.1", "GL000191.1", "GL000227.1", "GL000228.1", "GL000214.1", "GL000221.1",
                                 "GL000209.1", "GL000218.1", "GL000220.1", "GL000213.1", "GL000211.1", "GL000199.1",
                                 "GL000217.1", "GL000216.1", "GL000215.1", "GL000205.1", "GL000219.1", "GL000224.1",
                                 "GL000223.1", "GL000195.1", "GL000212.1", "GL000222.1", "GL000200.1", "GL000193.1",
                                 "GL000194.1", "GL000225.1", "GL000192.1", "NC_007605"}
        alignment_file = pysam.VariantFile(self.valid_vcf_file)
        obtained_sequence_names = self.vb_scanner.get_variant_sequence_names(alignment_file)
        alignment_file.close()
        self.assertSetEqual(obtained_sequence_names, sequence_names_answer, "The returned set of sequence names should "
                            f"have been {sequence_names_answer}")
