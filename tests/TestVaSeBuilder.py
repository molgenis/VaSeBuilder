#Import required libraries/modules
import unittest
import os

#Import required class
from VaSeBuilder import VaSeBuilder


#Unittest class for the VaSeBuilder class.
class TestVaSeBuilder(unittest.TestCase):
	
	#Set up requirements for this test class.
	def setUp(self):
		self.vsBuilder = VaSeBuilder()
	
	
	#Tests that two variant reads are obtained
	def test_getVariantReads_pos(self):
		answerList = ['SRR1039513.12406160', 'SRR1039513.12406160']
		varReads = vsBuilder.getVariantReads("16", 247990, "testdata/bamDir/vaseutest_1.bam")
		varReadNames = [x.query_name for x in varReads]
		assertEqual(varReadNames, answerList, "The lists should be identical but are not")
	
	#Tests that no reads are obtained
	def test_getVariantReads_neg(self):
		assertNone(self.vsBuilder.getVariantReads("16", 1, "testdata/bamDir/vaseutest_1.bam"))
	
	
	
	#Tests that the context of BAM reads associated to a variant is determined correctly.
	def test_determineContext_pos(self):
		
	
	def test_determineContext_neg(self):
		
	
	
	#Clean up.
	def tearDown(self):
		self.vsBuilder.dispose()