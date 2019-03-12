import unittest
from VariantContext import VariantContext


class TestVariantContext(unittest.TestCase):
	
	def setUp(self):
		self.variantContext = VariantContext('SNP21_1001000', '21', 100, 1000)
		self.donorReads = ['SRR21_001', 'SRR21_002', 'SRR21_003']
		self.acceptorReads = ['SRR21_100', 'SRR21_200', 'SRR21_300', 'SRR21_400']
	
	
	# Tests that the context variant id is indeed 'SNP21_1001000'
	def test_getVariantId(self):
		self.assertEqual(self.variantContext.getVariantId(), 'SNP21_1001000', 'The two variant identiers should be equal')
	
	
	# Tests that the context chromosome is indeed '21'
	def test_getContextChrom(self):
		self.assertEqual(self.variantContext.getContextChrom(), '21', 'The two chromosome names should have been 21')
	
	
	# Tests that the context start is indeed 100
	def test_getContextStart(self):
		self.assertEqual(self.variantContext.getContextStart(), 100, 'The two context start positions should have been 100')
	
	
	# Tests that the context end is indeed 1000
	def test_getContextEnd(self):
		self.assertEqual(self.variantContext.getContextEnd(), 1000, 'The two context end positions should have been 1000')
	
	
	# Tests that the context length is indeed 900
	def test_getContextLength(self):
		self.assertEqual(self.variantContext.getContextLength(), 900, 'The two contexts should have been 900 in length')
	
	
	
	# Tests that a single donor read identifier is added to the variant context
	def test_addDonorReadIdentifier(self):
		tmpContext = VariantContext('SNP21_1001000', '21', 100, 1000)
		tmpContext.addDonorReadIdentifier('tmp21_000')
		self.assertListEqual(tmpContext.contextDonorReadIds, ['tmp21_000'], 'The donor read identifier SRR21_000 should have been added')
	
	
	# Tests that the single acceptor read identifier is added to the variant context
	def test_addAcceptorReadIdentifier(self):
		tmpContext = VariantContext('SNP21_1001000', '21', 100, 1000)
		tmpContext.addAcceptorReadIdentifier('tmp21_000')
		self.assertListEqual(tmpContext.contextAcceptorReadIds, ['tmp21_000'], 'The acceptor read identifier SRR21_000 should have been added')
	
	
	# Tests that the list of donor read identifiers are indeed added to the variant context
	def test_addDonorReadIdentifiers(self):
		self.variantContext.addDonorReadIdentifiers(self.donorReads)
		self.assertListEqual(self.variantContext.contextDonorReadIds, self.donorReads, 'The list of read identifiers should have been saved in the context')
	
	
	# Tests that the list of acceptor read identifiers are indeed added to the variant context
	def test_addAcceptorReadIdentifiers(self):
		self.variantContext.addAcceptorReadIdentifiers(self.acceptorReads)
		self.assertListEqual(self.variantContext.contextAcceptorReadIds, self.acceptorReads, 'The list of acceptor read identifiers should have been saved in the context')
	
	
	# Tests that the number of donor reads in the variant context is indeed 3.
	def test_getNumberOfDonorReads(self):
		self.variantContext.addDonorReadIdentifiers(self.donorReads)
		self.assertEqual(self.variantContext.getNumberOfDonorReads(), 3, 'Three donor reads should have been in the variant context')
	
	
	# Tetst that the number of acceptor reads in the variant context is indeed 4.
	def test_getNumberOfAcceptorReads(self):
		self.variantContext.addAcceptorReadIdentifiers(self.acceptorReads)
		self.assertEqual(self.variantContext.getNumberOfAcceptorReads(), 4, 'Four acceptor reads should have been in the variant context')
	
	
	# Tests that the saved list of donor reads is the same as the one (self.donorReads) that we added
	def _test_getContextDonorReadIds(self):
		self.assertListEqual(self.variantContext.getContextDonorReadIds(), self.donorReads, 'Both list of donor read identifiers should have been the same')
	
	
	# Tests that the saved list of accepor reads is the same as the one (self.acceptorReads) that we added
	def _test_getContextAcceptorReadIds(self):
		self.assertListEqual(self.variantContext.getContextAcceptorReadIds(), self.acceptorReads, 'Both list of donor read identifiers should have been the same')
	
	
	# Tests that the read idenfier 'SRR21_001' is indeed added to the list of donor read identifiers
	def test_donorReadIsInVariantContext_pos(self):
		self.variantContext.addDonorReadIdentifiers(self.donorReads)
		self.assertTrue(self.variantContext.donorReadIsInVariantContext('SRR21_001'), 'Read identifier SRR21_001 should have been found')
	
	
	# Tests that the read identifier 'SRR21_000' is indeed not in the list of donor read identifiers
	def test_donorReadIsInVariantContext_neg(self):
		self.variantContext.addAcceptorReadIdentifiers(self.donorReads)
		self.assertFalse(self.variantContext.donorReadIsInVariantContext('SRR21_000'), 'Read identifier SRR21_000 should not have been found')
	
	
	# Tests that the read idenfier 'SRR21_100' is indeed added to the list of acceptor read identifiers
	def test_acceptorReadIsInVariantContext_pos(self):
		self.variantContext.addAcceptorReadIdentifiers(self.acceptorReads)
		self.assertTrue(self.variantContext.acceptorReadIsInVariantContext('SRR21_100'), 'Read identifier SRR21_100 should have been found')
	
	
	# Tests that the read identifier 'SRR21_000' is indeed not in the list of acceptor read identifiers
	def test_acceptorReadIsInVariantContext_neg(self):
		self.variantContext.addAcceptorReadIdentifiers(self.acceptorReads)
		self.assertFalse(self.variantContext.acceptorReadIsInVariantContext('SRR21_000'), 'Read identifier SRR21_000 should not have been found')


TestVariantContext()