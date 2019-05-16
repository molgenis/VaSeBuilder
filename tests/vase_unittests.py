#!/usr/bin/env python
# Import the unittest module
import sys
import unittest
sys.path.append('H:\Data\Repositories\VaSeBuilder')	# Temporarily append the path to the unittests for testing purposes.

# Import the unittest classes.
#import TestParamChecker
#import TestVcfBamScanner
import TestDonorBamRead
#import TestVaSeBuilder



# Create the the test loader and suite.
loader = unittest.TestLoader()
suite = unittest.TestSuite()


# Add the three test classes that test the functionality of each program class.
#suite.addTests(loader.loadTestsFromModule(TestParamChecker))
#suite.addTests(loader.loadTestsFromModule(TestVcfBamScanner))
suite.addTests(loader.loadTestsFromModule(TestDonorBamRead))
#suite.addTests(loader.loadTestsFromModule(TestVaSeBuilder))


# Run the tests.
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)