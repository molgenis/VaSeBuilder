import unittest

# Import the required VaSe classes
from DonorBamRead import DonorBamRead
from OverlapContext import OverlapContext
from VariantContext import VariantContext
from VariantContextFile import VariantContextFile

class TestVariantContextFile(unittest.TestCase):
	def setUp(self):
		
		# Create the variables containing context answers and values
		self.contextIdAnswer = '21_9411259'
		self.acceptorContextAnswer = OverlapContext()
		self.donorContextAnswer = OverlapContext()
		
		# Create the variables containing all the VariantContext answers
		self.variantContextAnswer = VariantContext()
		self.variantContextForSetOps = VariantContext()
		self.allVarconAcceptorReadsAnswer = []
		self.allVarconDonorReadsAnswer = []
		
		# Create the variables containg the VariantContextFile answers
		self.variantContextsAnswer = {self.contextIdAnswer:self.variantContextAnswer}
		self.posFilterValToUse = 'aap'
		self.negFilterValToUse = 'jan'
		self.filterListToUse = ['aap', 'noot', 'mies']
		
		# Create the variables containing the answer of the VariantContextFile
		self.variantContextFile = VariantContextFile()
	
	
	
	# ====================PERFORMS THE TESTS FOR THE GETTER METHODS====================
	def test_getVariantContexts(self):
		self.assertDictEqual(self.variantContextFile.getVariantContexts(), self.variantContextsAnswer, "The returned variant contexts are not what was expected")
	
	def test_getVariantContext(self):
		self.assertEqual(self.variantContextFile.getVariantContext(self.contextIdAnswer).toString(), self.variantContextAnswer.toString(), "The returned variant context is not what was expected")
	
	def test_getVariantContext_None(self):
		self.assertIsNone(self.variantContextFile.getVariantContext('22_9411255'), "The requested variant context should not have existed and should have therefore been None")
	
	def test_getAcceptorContext(self):
		self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer).toString(), self.acceptorContextAnswer.toString(), "The returned acceptor contexts is not what was expected")
	
	def test_getAcceptorContext_None(self):
		self.assertIsNone(self.variantContextFile.getAcceptorContext('22_9411255'), "The requested acceptor context should not have existed and should have therefore been None")
	
	def test_getDonorContext(self):
		self.assertEqual(self.variantContextFile.getAcceptorContext(self.contextIdAnswer), self.donorContextAnswer.toString(), "The returned donor context is not what was expected")
	
	def test_getDonorContext_None(self):
		self.assertIsNone(self.variantContextFile.getDonorContext('22_9411255'), "The requested donor context should not have existed and should have therefore been None")
	
	def test_getAllVariantContextAcceptorReads(self):
		self.assertListEqual(self.variantContextFile.getAllVariantContextAcceptorReads(), self.allVarconAcceptorReadsAnswer, "")
	
	def test_getAllVariantContextDonorReads(self):
		self.assertListEqual(self.variantContextFile.getAllVariantContextAcceptorReads(), self.allVarconDonorReadsAnswer, "")
	
	def test_getAllVariantContextAcceptorReadIds(self):
	
	
	def test_getAllVariantContextDonorReadIds(self):
		
	
	
	# ====================PERFORM THE TESTS FOR READING A VARIANT CONTEXT FILE====================
	#def test_readVariantContextFile_pos(self)
	
	def test_readVariantContextFile_neg(self):
		with self.assertRaises(IOError) as cm:
			self.variantContextFile.readVariantContextFile('nonexisting.txt')
		raisedException = cm.exception
		self.assertEqual(raisedException.error_code, 3)
	
	def test_passesFilter_pos(self):
		self.assertTrue(self.variantContextFile.passesFilter(self.posFilterValToUse), f"The value {self.posFilterValToUse} should have been in the filter list and therefore return True")
		
	def test_passesFilter_neg(self):
		self.assertFalse(self.variantContextFile.passesFilter(self.negFilterValToUse), f"The value {self.negFilterValToUse} should not have been in the filter list and therefore return False")
	
	
	
	# ====================PERFORM THE TESTS FOR SET OPERATIONS ON TWO VARIANT CONTEXT FILES====================
	def getVariantContextsUnion_pos(self, otherVarconFile)
	def getVariantContextsUnion_neg(self, otherVarconFile)
	def getVariantContextsIntersect_pos(self, otherVarconFile)
	def getVariantContextsIntersect_neg(self, otherVarconFile)
	def getVariantContextsDifference_pos(self, otherVarconFile)
	def getVariantContextsDifference_neg(self, otherVarconFile)
	def getVariantContextsSymmetricDifference_pos(self, otherVarconFile)
	def getVariantContextsSymmetricDifference_neg(self, otherVarconFile)
