#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os
import sys
import pysam

# Import required class
sys.path.insert(0,'./../')
from VcfBamScanner import VcfBamScanner


# Unittest class for the VcfBamScanner class
class TestVcfBamScanner(unittest.TestCase):
    
    # Set up some things for the testing methods.
    def setUp(self):
        self.vb_scanner = VcfBamScanner()
        self.vcf_folders = ['testdata/vcfDir']
        self.bam_folders = ['testdata/bamDir']
        self.none_folders = ['testdata/doesnotexist']
    
    
    
    # Tests the scanning of VCF folders that exist and contain vcf/vcf.gz files
    def test_scanVcfFolders_pos(self):
        answer_dict = {'SRR1039508': 'testdata/vcfDir/SRR1039508.vcf', 'SRR1039512': 'testdata/vcfDir/SRR1039512.vcf'}
        result_dict = self.vb_scanner.scan_vcf_folders(self.vcf_folders)
        self.assertDictEqual(result_dict, answer_dict, "Dicts should have been the same")
    
    
    
    # Tests the scanning of VCF folders that do not exist / contain no vcf/vcf.gz files.
    def test_scanVcfFolders_neg(self):
        answer_dict = {}
        result_dict = self.vb_scanner.scan_vcf_folders(self.none_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been empty but are not")
    
    
    
    # Tests the scanning of folders containg BAM files.
    def test_scanBamFolders_pos(self):
        answer_dict = {'SRR1039508': 'testdata/bamDir/SRR1039508.bam', 'SRR1039512': 'testdata/bamDir/SRR1039512.bam'}
        result_dict = self.vb_scanner.scan_bam_folders(self.bam_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been the same")
    
    
    
    # Tests the scanning of folders that do not exist / contain no BAM files.
    def test_scanBamFolders_neg(self):
        answer_dict = {}
        result_dict = self.vb_scanner.scan_bam_folders(self.none_folders)
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been empty")
    
    
    
    # Tests that a BAM file has sample information.
    def test_bamHasSampleName_pos(self):
        bam_file = pysam.AlignmentFile("testdata/bamDir/SRR1039508.bam")
        result_bool = self.vb_scanner.bam_has_sample_name(bam_file)
        bam_file.close()
        self.assertTrue(result_bool)
    
    
    
    # Tests that a BAM file has no sample information.
    def test_bamHasSampleName_neg(self):
        bam_file = pysam.AlignmentFile("testdata/noSampleDir/noSampleBam.bam")
        result_bool = self.vb_scanner.bam_has_sample_name(bam_file)
        bam_file.close()
        self.assertFalse(result_bool)
    
    
    # Tests whether the VCF to BAM map will be constructed properly.
    def test_getVcfToBamMap(self):
        answer_dict = {'testdata/vcfDir/SRR1039508.vcf': 'testdata/bamDir/SRR1039508.bam', 'testdata/vcfDir/SRR1039512.vcf': 'testdata/bamDir/SRR1039512.bam'}
        vcf_fs = self.vb_scanner.scan_vcf_folders(self.vcf_folders)    # Provide the VcfBamScanner with a valid list of VCF files
        bam_fs = self.vb_scanner.scan_bam_folders(self.bam_folders)    # Provide the VcfBamScanner with a valid list of BAM files.
        result_dict = self.vb_scanner.get_vcf_to_bam_map()
        self.assertDictEqual(result_dict, answer_dict, "Both dicts should have been the same")
