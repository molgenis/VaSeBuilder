#!'/usr/bin/env python
import unittest
from VariantContextFile import VariantContextFile
from VariantContext import VariantContext

class TestVariantContextFile(unittest.TestCase):
	
	# Create all the variables that will be used in testing
	def setUp(self):
		self.testNonExistingContextId = 'SNP21_000'
		self.testValidVarConId = 'SNP21_0100'
		self.testValidVarConChr = '21'
		self.testValidVarConStart = 0
		self.testValidVarConEnd = 100
		self.testVariantContextObj = VariantContext(self.testValidVarConId, self.testValidVarConChr, self.testValidVarConStart, self.testValidVarConEnd)
		self.testVariantContexts = {'SNP21_0100' : VariantContext('SNP21_0100', '21', 0, 100),
			'SNP21_1001000' : VariantContext('SNP21_1001000', '21', 100, 1000)}
		self.variantContextFile = VariantContextFile('test_varcon.txt')
	
	
	# Tests that the file is read properly
	def test_readVariantContextFile(self):
		#self.variantContextFile = VariantContextFile('test_varcon.txt')
		self.assertDictEqual(self.variantContextFile.variantContexts, self.testVariantContexts)
	
	# Tests that the variant context have been read
	def test_getVariantContexts(self):
		self.assertDictEqual(self.variantContextFile.getVariantContexts(), self.testVariantContexts, '')
	
	# Tests that a variant context is indeed in the variant context file
	def test_hasVariantContext_pos(self):
		self.assertTrue(self.variantContextFile.hasVariantContext(self.testValidVarConId), "Variant ID "+str(self.testValidVarConId)+" should have been in the variant context file")
	
	# Tests that a variant context is indeed 
	def test_hasVariantContext_neg(self):
		self.assertFalse(self.variantContextFile.hasVariantContext(self.testNonExistingContextId), "Variant ID "+str(self.testValidVarConId)+" should not have been in the variant context file")
	
	# Tests that an existing variant context is returned.
	def test_getVariantContext_pos(self):
		varContext = self.variantContextFile.getVariantContext(self.testValidVarConId)
		self.assertEqual(varContext.getVariantId(), self.testVariantContextObj.getVariantId(), "The variant context IDs should both have been "+str(self.testValidVarConId))
		self.assertEqual(varContext.get_context_chrom(), self.testVariantContextObj.getContextChrom(), "The variant context chroms should both have been " + str(self.testValidVarConChr))
		self.assertEqual(varContext.get_context_start(), self.testVariantContextObj.getContextStart(), "The variant context start positions should both have been " + str(self.testValidVarConStart))
		self.assertEqual(varContext.get_context_end(), self.testVariantContextObj.getContextEnd(), "The variant context end positions should both have been " + str(testValidVarConEnd))
	
	# Tests that a non existent variant context is not returned
	def test_getVariantContext_neg(self):
		self.assertIsNone(self.variantContextFile.getVariantContext(self.testNonExistingContextId), "Searching for variant context "+str(self.testNonExistingContextId)+ " should not have been found")
	
	# Tests that the number of variant contexts is indeed 2
	def test_getNumberOfContexts(self):
		self.assertEqual(self.variantContextFile.getNumberOfContexts(), 2, 'The number of variant contexts should both have been 2')
	
