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
		self.vbScanner = VcfBamScanner()
		self.vcfFolders = ['testdata/vcfDir']
		self.bamFolders = ['testdata/bamDir']
		self.noneFolders = ['testdata/doesnotexist']
	
	
	
	# Tests the scanning of VCF folders that exist and contain vcf/vcf.gz files
	def test_scanVcfFolders_pos(self):
		answerDict = {'SRR1039508': 'testdata/vcfDir/SRR1039508.vcf', 'SRR1039512': 'testdata/vcfDir/SRR1039512.vcf'}
		resultDict = self.vbScanner.scanVcfFolders(self.vcfFolders)
		self.assertDictEqual(resultDict, answerDict, "Dicts should have been the same")
	
	
	
	# Tests the scanning of VCF folders that do not exist / contain no vcf/vcf.gz files.
	def test_scanVcfFolders_neg(self):
		answerDict = {}
		resultDict = self.vbScanner.scanVcfFolders(self.noneFolders)
		self.assertDictEqual(resultDict, answerDict, "Both dicts should have been empty but are not")
	
	
	
	# Tests the scanning of folders containg BAM files.
	def test_scanBamFolders_pos(self):
		answerDict = {'SRR1039508': 'testdata/bamDir/SRR1039508.bam', 'SRR1039512': 'testdata/bamDir/SRR1039512.bam'}
		resultDict = self.vbScanner.scanBamFolders(self.bamFolders)
		self.assertDictEqual(resultDict, answerDict, "Both dicts should have been the same")
	
	
	
	# Tests the scanning of folders that do not exist / contain no BAM files.
	def test_scanBamFolders_neg(self):
		answerDict = {}
		resultDict = self.vbScanner.scanBamFolders(self.noneFolders)
		self.assertDictEqual(resultDict, answerDict, "Both dicts should have been empty")
	
	
	
	# Tests that a BAM file has sample information.
	def test_bamHasSampleName_pos(self):
		bamFile = pysam.AlignmentFile("testdata/bamDir/SRR1039508.bam")
		resultBool = self.vbScanner.bamHasSampleName(bamFile)
		bamFile.close()
		self.assertTrue(resultBool)
	
	
	
	# Tests that a BAM file has no sample information.
	def test_bamHasSampleName_neg(self):
		bamFile = pysam.AlignmentFile("testdata/noSampleDir/noSampleBam.bam")
		resultBool = self.vbScanner.bamHasSampleName(bamFile)
		bamFile.close()
		self.assertFalse(resultBool)
	
	
	# Tests whether the VCF to BAM map will be constructed properly.
	def test_getVcfToBamMap(self):
		answerDict = {'testdata/vcfDir/SRR1039508.vcf': 'testdata/bamDir/SRR1039508.bam', 'testdata/vcfDir/SRR1039512.vcf': 'testdata/bamDir/SRR1039512.bam'}
		vcfFs = self.vbScanner.scanVcfFolders(self.vcfFolders)	# Provide the VcfBamScanner with a valid list of VCF files
		bamFs = self.vbScanner.scanBamFolders(self.bamFolders)	# Provide the VcfBamScanner with a valid list of BAM files.
		resultDict = self.vbScanner.getVcfToBamMap()
		self.assertDictEqual(resultDict, answerDict, "Both dicts should have been the same")
