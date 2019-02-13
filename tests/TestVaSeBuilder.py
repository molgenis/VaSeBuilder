#!/usr/bin/env python
# Import required libraries/modules
import unittest
from unittest import mock
from unittest.mock import patch, mock_open, call
import os
import sys
import pysam

# Import required class
sys.path.insert(0, './../')
from VaSeBuilder import VaSeBuilder


# Unittest class for the VaSeBuilder class.
class TestVaSeBuilder(unittest.TestCase):
	
	# Set up requirements for this test class.
	def setUp(self):
		self.vsBuilder = VaSeBuilder("aap")
		self.varConMap = {'SNP16_247990':['16', 247986, 56508478]}
	
	
	# Tests that two variant reads are obtained
	def test_getVariantReads_pos(self):
		answerList = ['SRR1039513.12406160', 'SRR1039513.12406160']
		varReads = self.vsBuilder.getVariantReads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
		varReadNames = [x.query_name for x in varReads]
		self.assertListEqual(varReadNames, answerList, "The lists should be identical but are not")
	
	# Tests that no reads are obtained
	def test_getVariantReads_neg(self):
		self.assertListEqual(self.vsBuilder.getVariantReads("16", 1, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb')), [], "Both list should be empty")
	
	
	
	# Tests that the context of BAM reads associated to a variant is determined correctly.
	def test_determineContext_pos(self):
		answerList = ["16", 247985, 56508477]
		varReads = self.vsBuilder.getVariantReads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
		resultList = self.vsBuilder.determineContext(varReads)
		self.assertListEqual(resultList, answerList, "The contexts should have been the same.")
		
	# Test that the context of none existing BAM reads has no context.
	def test_determineContext_neg(self):
		answerList = []
		varReads = self.vsBuilder.getVariantReads("16", 1, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
		resultList = self.vsBuilder.determineContext(varReads)
		self.assertListEqual(resultList, answerList, "")
	
	
	
	# Tests that a variant at position 250000 is within the context [247986 .. 56508478]
	def test_isInContext_pos(self):
		self.vsBuilder.variantContextMap = self.varConMap
		self.assertTrue(self.vsBuilder.isInContext('16', 250000))
		self.vsBuilder.variantContextMap = {}	# Reset the variantContextMap of the vsBuilder object.
	
	# Tests that a variant at position 247985 is not in the context [247986 .. 56508478]
	def test_isInContext_neg(self):
		self.vsBuilder.variantContextMap = self.varConMap
		self.assertFalse(self.vsBuilder.isInContext('16', 247985))
		self.vsBuilder.variantContextMap = {}	# Reset the variantContextMap of the vsBuilder object.
	
	
	
	#Tests that a variant has already been processed.
	def test_variantAlreadyProcessed_pos(self):
		self.vsBuilder.variantContextMap = {}
		self.vsBuilder.variantContextMap = self.varConMap
		self.assertTrue(self.vsBuilder.variantAlreadyProcessed('SNP16_247990'))
	
	#Tests that a variant hasn't been processed.
	def test_variantAlreadyProcessed_neg(self):
		self.vsBuilder.variantContextMap = {}
		self.vsBuilder.variantContextMap = self.varConMap
		self.assertFalse(self.vsBuilder.variantAlreadyProcessed('SNP16_250000'))
	
	
	
	#Tests that a read is the required read.
	def test_isRequiredRead_pos(self):
		varReads = self.vsBuilder.getVariantReads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
		if(varReads[1].is_read1):
			self.assertTrue(self.vsBuilder.isRequiredRead(varReads[1], 'F'), "Should have been true for forward read")
		else:
			self.assertTrue(self.vsBuilder.isRequiredRead(varReads[0], 'F'), "Should have been true for forward read")
	
	#Test that a read is not the required read.
	def test_isRequiredRead_neg(self):
		varReads = self.vsBuilder.getVariantReads("16", 247990, pysam.AlignmentFile("testdata/valbam/SRR1039513.bam", 'rb'))
		if(varReads[0].is_read2):
			self.assertFalse(self.vsBuilder.isRequiredRead(varReads[0], "F"), "Should have been false for forward read")
		else:
			self.assertFalse(self.vsBuilder.isRequiredRead(varReads[1], "F"), "Should have been false for forward read")
	
	
	
	#Tests that a saved context for a specified variant is indeed returned.
	def test_getVariantContext_pos(self):
		self.vsBuilder.variantContextMap = self.varConMap
		resultList = self.vsBuilder.getVariantContext('SNP16_247990')
		self.vsBuilder.variantContextMap = {}
		self.assertListEqual(resultList, ['16', 247986, 56508478], "Contexts should have been equal")
	
	#Tests that a None context is returned for a non existent variant.
	def test_getVariantContext_neg(self):
		self.vsBuilder.variantContextMap = self.varConMap
		resultList = self.vsBuilder.getVariantContext('SNP15_10000')
		self.vsBuilder.variantContextMap = {}
		self.assertNone(resultList)
	
	
	
	#Tests that specific datais indeed written to a file.
	def test_writeVariantsContexts(self):
		varconMapToWrite = {'SNP16_001':['16', 10, 100], 'SNP16_002':['16', 20, 200]}
		resultCalls =[call("Variant\tChrom\tStart\tStop\n"), call("SNP16_001\t16\t10\t100\n"), call("SNP16_002\t16\t20\t200\n")]
		m = mock.mock_open()
		with mock.patch(self.vsBuilder.writeVariantsContexts.open , m, create=True) as mocked_open:
			self.vsBuilder.writeVariantsContexts()
			file_handle = mocked_open()
			file_handle.write.assert_has_calls(resultCalls)
	
	
	#Tests that a fastQ is build.
	def test_buildFastQ_pos(self):
		print("Implementing...")
	
	#Tests that a fastQ is not build.
	def test_buildFastQ_neg(self):
		print("Implementing...")
