#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os

# Import required class
from ParamChecker import ParamChecker


# Unittest class for the ParamChecker test.
class TestParamChecker(unittest.TestCase):
	
	# Set up required 
	def setUp(self):
		self.paramCheck = ParamChecker()
		self.paramList = {'vcfin': 'testdata/vcfDir',
			'bamin': 'testdata/bamDir',
			'valbam': 'testdata/valbam.bam',
			'valfastq1': 'testdata/fqDir/vaseval_F.fq.gz',
			'valfastq2': 'testdata/fqDir/vaseval_R.fq.gz',
			'fastqout': 'testdata/outDir',
			'varcon': 'testdata/outDir/varcon.txt',
			'varbread': 'testdata/outDir/varbread.txt',
			'nistbread': 'testdata/outDir/nistbread.txt'
		}
	
	
	# Tests the 
	def test_checkLog_pos(self):
		
	
	# Tests the 
	def test_checkLog_neg(self):
		
	
	
	
	def test_checkFolderContents_pos(self):
		
	
	def test_checkFolderContents_neg(self):
		
	
	
	
	def test_checkFoldersExist_pos(self):
		
	
	
	def test_checkFoldersExist_neg(self):
		
	
	
	
	def test_checkFileExists_pos(self):
	def test_checkFileExists_neg(self):
	
	
	
	
	# Test that the output location for out.txt is valid.
	def test_isValidOutputLocation_pos(self):
		assertTrue(self.paramCheck.isValidOutputLocation("testdata/outDir/out.txt"))
	
	# Test that the output location for out.txt is invalid.
	def test_isValidOutputLocation_neg(self):
		assertFalse(self.paramCheck.isValidOutputLocation("testdata/doesnotexist/out.txt"))
	
	
	
	# Tests that the foldername a file is located in gets returned properly.
	def test_getFolderName(self):
		assertEqual( self.paramCheck.getFolderName("testdata/noSampleDir/noSampleBam.bam"), "noSampleDir", "Folder names should be equal bur are not")
	
	
	
	# Tests that all provided parameters are ok. First test should return True, all others false due to one parameter missing
	def test_checkParameters_pos(self):
		assertTrue(self.paramCheck(self.paramList))
	
	def test_checkParameters_noVcfIn(self):
		parList = self.paramList.copy()
		parList['vcfin'] = 'testdata/doesnotexist'	# Set the list VCF folders to a folder that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noBamIn(self):
		parList = self.paramList.copy()
		parList['bamin'] = 'testdata/doesnotexist'	# Set the list BAM folders to a folder that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noValBam(self):
		parList = self.paramList.copy()
		parList['valbam'] = 'testdata/doesnotexist.bam'	# Set the location of the bam file to one that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noValFq1(self):
		parList = self.paramList.copy()
		parList['valfastq1'] = 'testdata/doesnotexist.fq'	# Set the location of the fastq file to one that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noValFq2(self):
		parList = self.paramList.copy()
		parList['valfastq2'] = 'testdata/doesnotexist.fq'	# Set the location of the fastq file to one that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noFqOut(self):
		parList = self.paramList.copy()
		parList['fastqout'] = 'testdata/doesnotexist/faqOut.fq'	# Set the output location of the fastq file to a folder that does not exist
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noVarcon(self):
		parList = self.paramList.copy()
		parList['varcon'] = 'testdata/doesnotexist/varcon.txt'	# Set the location for the varcon output file to a folder that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noVarBread(self):
		parList = self.paramList.copy()
		parList['varbread'] = 'testdata/doesnotexist/varbread.txt'	# Set the location for the varbread output file to a folder that does not exist.
		assertFalse(self.paramCheck(parList))
	
	def test_checkParameters_noNistBread(self):
		parList = self.paramList.copy()
		parList['nistbread'] = 'testdata/doesnotexist/nistbread.txt'	# Set the location for the nistbread output file to a folder that does not exist.
		assertFalse(self.paramCheck(parList))
	
	
	
	# Clean up
	def tearDown(self):
		self.paramCheck.dispose()
		self.paramList.dispose()
