#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os
import pysam

# Import required class
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
		answerDict = {'VaSeUTest1': 'testdata/vcfDir/vaseutest_1.vcf.gz', 'VaSeUTest2': 'testdata/vcfDir/vaseutest_2.vcf.gz'}
		resultDict = self.vbScanner.scanVcfFolders(self.vcfFolders)
		assertDictEqual(resultDict, answerDict, "Dicts should have been the same")
	
	
	
	# Tests the scanning of VCF folders that do not exist / contain no vcf/vcf.gz files.
	def test_scanVcfFolders_neg(self):
		answerDict = {}
		resultDict = self.vbScanner.scanVcfFolders(self.noneFolders)
		assertDictEqual(resultDict, answerDict, "Both dicts should have been empty but are not")
	
	
	
	# Tests the scanning of folders containg BAM files.
	def test_scanBamFolders_pos(self):
		answerDict = {'VaSeUTest1': 'testdata/bamDir/vaseutest_1.bam', 'VaSeUTest2': 'testdata/bamDir/vaseutest_2.bam'}
		resultDict = self.vbScanner.scanBamFolders(self.bamFolders)
		assertDictEqual(resultDict, answerDict, "Both dicts should have been the same")
	
	
	
	# Tests the scanning of folders that do not exist / contain no BAM files.
	def test_scanBamFolders_neg(self):
		answerDict = {}
		resultDict = self.vbScanner.scanBamFolders(self.noneFolders)
		assertDictEqual(resultDict, answerDict, "Both dicts should have been empty")
	
	
	
	# Tests that a BAM file has sample information.
	def test_bamHasSampleName_pos(self):
		bamFile = pysam.AlignmentFile("testdata/bamDir/vasesutest_1.bam")
		resultBool = vbScanner.bamHasSampleName(bamFile)
		bamFile.close()
		assertFalse(resultBool)
	
	
	
	# Tests that a BAM file has no sample information.
	def test_bamHasSampleName_neg(self):
		bamFile = pysam.AlignmentFile("testdata/noSampleDir/noSampleBam.bam")
		resultBool = vbScanner.bamHasSampleName(bamFile)
		bamFile.close()
		assertFalse(resultBool)
	
	
	# Tests whether the VCF to BAM map will be constructed properly.
	def test_getVcfToBamMap(self):
		answerDict = {'testdata/vcfDir/vaseutest_1.vcf.gz': 'testdata/bamDir/vaseutest_1.bam', 'testdata/vcfDir/vaseutest_1.vcf.gz': 'testdata/bamDir/vaseutest_2.bam'}
		vcfFs = vbScanner.scanVcfFolders(self.vcfFolders)	# Provide the VcfBamScanner with a valid list of VCF files
		bamFs = vbScanner.scanBamFolders(self.bamFolders)	# Provide the VcfBamScanner with a valid list of BAM files.
		resultDict = vbScanner.getVcfToBamMap()
		assertDictEqual(resultDict, answerDict, "Both dicts should have been the same")
	
	
	
	# Clean up variables and such.
	def tearDown(self):
		self.vbScanner.dispose()
		self.vcfFolders.dispose()
		self.bamFolders.dispose()
		self.noneFolders.dispose()
