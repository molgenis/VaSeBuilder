#!/usr/bin/env python
import unittest
from GvcfVariant import GvcfVariant

class TestGvcfVariant(unittest.TestCase):
	
	# Use the following data as test:
	#21	9411327	.	C	G	257.77	.	AC=1;AF=0.5;AN=2;BaseQRankSum=0.222;ClippingRankSum=-0.37;DP=52;ExcessHet=3.0103;FS=0;MLEAC=1;MLEAF=0.5;MQ=54.26;MQRankSum=-0.271;QD=9.55;ReadPosRankSum=0.025;SOR=1.179	GT:AD:DP:GQ:PL	0/1:16,11:27:99:286,0,438
	def setUp(self):
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
		
		self.gvcfVariant = GvcfVariant(self.varId, self.varChrom, self.varPos, self.varRef, self.varAlt, self.varQual, self.varFilter, self.varInfo, self.varFormatSample)
	
	
	# Tests that the variant identifier '.' is indeed returned
	def test_getVariantId(self):
		self.assertEqual(self.gvcfVariant.getVariantId(), self.varId, "The two variant identifiers should have been a \'.\'")
	
	
	def test_getChromPosId(self):
		self.assertEqual(self.gvcVariant.getChromPosId(), self.chromPosId, 'The two chrompos identifiers should have been SNP21_9411327')