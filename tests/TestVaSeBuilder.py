#!/usr/bin/env python
# Import required libraries/modules
import unittest
import os

# Import required class
from VaSeBuilder import VaSeBuilder


# Unittest class for the VaSeBuilder class.
class TestVaSeBuilder(unittest.TestCase):
	
	# Set up requirements for this test class.
	def setUp(self):
		self.vsBuilder = VaSeBuilder()
	
	
	# Tests that two variant reads are obtained
	def test_getVariantReads_pos(self):
		answerList = ['SRR1039513.12406160', 'SRR1039513.12406160']
		varReads = vsBuilder.getVariantReads("16", 247990, "testdata/bamDir/vaseutest_1.bam")
		varReadNames = [x.query_name for x in varReads]
		assertEqual(varReadNames, answerList, "The lists should be identical but are not")
	
	# Tests that no reads are obtained
	def test_getVariantReads_neg(self):
		assertNone(self.vsBuilder.getVariantReads("16", 1, "testdata/bamDir/vaseutest_1.bam"))
	
	
	
	# Tests that the context of BAM reads associated to a variant is determined correctly.
	def test_determineContext_pos(self):
		answerList = ["16", 247986, 56508478]
		varReads = self.vsBuilder.getVariantReads("16", 247990, "testdata/bamDir/vaseutest_1.bam")
		resultList = self.vsBuilder.determineContext(varReads)
		assertListEqual(resultList, answerList, "The contexts should have been the same.")
		
	# Test that the context of none existing BAM reads has no context.
	def test_determineContext_neg(self):
		answerList = []
		varReads = self.vsBuilder.getVariantReads("16", 1, "testdata/bamDir/vaseutest_1.bam")
		resultList = self.vsBuilder.determineContext(varReads)
		assertListEqual(resultList, answerList, "")
	
	
	# Clean up.
	def tearDown(self):
		self.vsBuilder.dispose()