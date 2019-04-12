#!/usr/bin/env python
import unittest
from GvcfVariant import GvcfVariant

class TestGvcfVariant(unittest.TestCase):
	
	# Use the following data as test:
	#21	9411327	.	C	G	257.77	.	AC=1;AF=0.5;AN=2;BaseQRankSum=0.222;ClippingRankSum=-0.37;DP=52;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=54.26;MQRankSum=-0.271;QD=9.55;ReadPosRankSum=0.025;SOR=1.179	GT:AD:DP:GQ:PL	0/1:16,11:27:99:286,0,438
	def setUp(self):
		self.testFormatFieldPos = 'GT'
		self.testFieldNeg = 'NE'
		
		self.varId = '.'
		self.chromPosId = 'SNP21_9411327'
		self.varChrom = '21'
		self.varPos = 9411327
		self.varRef = 'G'
		self.varAlt = 'C'
		self.varQual = 257.77
		self.varFilter = '.'
		
		self.varInfo = {'AC' : 1.0,
			'AF' : 0.5,
			'AN' : 2.0,
			'BaseQRankSum' : 0.222,
			'ClippingRankSum' : -0.37,
			'DP': 52.0,
			'ExcessHet' : 3.0103,
			'FS' : 0.0,
			'MLEAC' : 1.0,
			'MLEAF' : 0.5,
			'MQ' : 54.26,
			'MQRankSum' : -0.271,
			'QD' : 9.55,
			'ReadPosRankSum' : 0.025,
			'SOR' : 1.179}
		
		self.varFormatSample = {'GT' : '0/1',
			'AD' : '16,11',
			'DP' : '27',
			'GQ' : '99',
			'PL' : '286,0,438'}
		
		self.calledInfo = {'gt' : '1/1',
			'ref' : 'G',
			'alt' : 'C',
			'pos' : '9411327',
			'res' : 'ci'}
		
		self.gvcfVariant = GvcfVariant(self.varId, self.varChrom, self.varPos, self.varRef, self.varAlt, self.varQual, self.varFilter, self.varInfo, self.varFormatSample)
	
	
	# Tests that the variant identifier '.' is indeed returned
	def test_getVariantId(self):
		self.assertEqual(self.gvcfVariant.getVariantId(), self.varId, "The two variant identifiers should have been "+str(self.varId))
	
	# Tests that the variant chrompos identifier is indeed 'SNP21_9411327'
	def test_getChromPosId(self):
		self.assertEqual(self.gvcfVariant.getChromPosId(), self.chromPosId, "The two chrompos identifiers should have been"+str(self.chromPosId))
	
	# Tests that the variant position is indeed 9411327
	def test_getVariantPos(self):
		self.assertEqual(self.gvcfVariant.getVariantPos(), self.varPos, "The two variant positions should have been "+str(self.varPos))
	
	# Tests that the variant reference is indeed a 'G'
	def test_getVariantRef(self):
		self.assertEqual(self.gvcfVariant.getVariantRef(), self.varRef, "The two variant references should have been a "+str(self.varRef))
	
	# Tests that the variant alternative is indeed a 'C'
	def test_getVariantAlt(self):
		self.assertEqual(self.gvcfVariant.getVariantAlt(), self.varAlt, "The two variant alternatives should have been a "+str(self.varAlt))
	
	# Tests that the variant qual is indeed 257.77
	def test_getVariantQual(self):
		self.assertEqual(self.gvcfVariant.getVariantQual(), self.varQual, "The two variant qualities should have been "+str(self.varQual))
	
	# Tests that the variant filter is set correctly.
	def test_getVariantFilter(self):
		self.assertEqual(self.gvcfVariant.getVariantFilter(), self.varFilter, "The two variant filters should have been "+str(self.varFilter))
	
	# Tests that the variant info is set correctly.
	def test_getVariantInfo(self):
		self.assertDictEqual(self.gvcfVariant.getVariantInfo(), self.varInfo, 'The two variant info dicts should have been the same')
	
	# Tests that the variant format and sample info has been set correctly
	def test_getVariantFormatSampleData(self):
		self.assertDictEqual(self.gvcfVariant.getVariantFormatSampleData(), self.varFormatSample, 'The two format-sample dicts should have been the same')
	
	
	# Tests that the variant info field 'AF' is indeed in the info dict
	def test_hasVariantInfoField_pos(self):
		self.assertTrue(self.gvcfVariant.hasVariantInfoField('AF'), 'The info field AF should have been in the variant info dict')
	
	# Tests that the variant info field 'NE' is indeed not in the info dict
	def test_hasVariantInfoField_neg(self):
		self.assertFalse(self.gvcfVariant.hasVariantInfoField(self.testFieldNeg), 'The info field NE should not have been in the variant info dict')
	
	
	# Tests that the variant has the format-sample field 'GT' in the format-sample dict
	def test_sampleHasFormatField_pos(self):
		self.assertTrue(self.gvcfVariant.sampleHasFormatField(self.testFormatFieldPos), 'The format-sample field GT should have been in the format-sample dict')
	
	# Tests that the variant does not have the format-sample field 'NE' in the format-sample dict
	def test_sampleHasFormatField_neg(self):
		self.assertFalse(self.gvcfVariant.sampleHasFormatField(self.testFieldNeg), 'The format-sample field NE should not have been in the format-sample dict')
	
	
	# Tests that sample genotype information is returned correctly
	def test_getSampleInfo_pos(self):
		self.assertEqual(self.gvcfVariant.getSampleInfo(self.testFormatFieldPos), '0/1', 'The returned sample genotype should have been 0/1')
	
	# Tests that a None value is returned for a non existing sample format value
	def test_getSampleInfo_neg(self):
		self.assertIsNone(self.gvcfVariant.getSampleInfo(self.testFieldNeg), 'No sample info should have been returned')
	
	
	# Tests that pipeline called variant info has not been set
	def test_hasCalledInfo_neg(self):
		self.assertFalse(self.gvcfVariant.hasCalledInfo(), 'There should have been no pipeline called variant info')
	
	# Tests that the called info is indeed set
	def test_setCalledInfo(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertDictEqual(self.gvcfVariant.calledInfo, self.calledInfo, 'The called info should have been set')
	
	# Tests that pipeline called variant info has been set
	def test_hasCalledInfo_pos(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertTrue(self.gvcfVariant.hasCalledInfo(), 'There should have been pipeline called variant info')
	
	# Tests that getCalledInfo returns the same dict as self.calledInfo
	def test_getCalledInfo(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertDictEqual(self.gvcfVariant.getCalledInfo(), self.calledInfo, 'The called info should have been the same')
	
	
	# Tests that the correct pipeline called genotype is returned
	def test_getCalledGt(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertEqual(self.gvcfVariant.getCalledGt(), '1/1', '')
	
	# Tests that the correct pipeline called reference allele is returned
	def test_getCalledRef(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertEqual(self.gvcfVariant.getCalledRef(), 'G', '')
	
	# Tests that the correct pipeline called alternative allele is returned
	def test_getCalledAlt(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertEqual(self.gvcfVariant.getCalledAlt(), 'C', '')
	
	# Tests that the correct pipeline called position is returned
	def test_getCalledPos(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertEqual(self.gvcfVariant.getCalledPos(), '9411327', '')
	
	# Tests that the correct pipeline called result is returned
	def test_getCalledResult(self):
		self.gvcfVariant.setCalledInfo(self.calledInfo)
		self.assertEqual(self.gvcfVariant.getCalledResult(), 'ci', '')
