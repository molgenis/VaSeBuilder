#!/usr/bin/env python
# Import the unittest module
import unittest

# Import the unittest classes.
import TestParamChecker
import TestVcfBamScanner
import TestVaSeBuilder


# Create the the test loader and suite.
loader = unittest.TestLoader()
suite = unittest.TestSuite()


# Add the three test classes that test the functionality of each program class.
suite.addTests(loader.loadTestsFromModule(TestParamChecker))
suite.addTests(loader.loadTestsFromModule(TestVcfBamScanner))
suite.addTests(loader.loadTestsFromModule(TestVaSeBuilder))


# Run the tests.
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)