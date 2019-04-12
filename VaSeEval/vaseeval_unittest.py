#!/usr/bin/env python
# Import the unittest module
import unittest

# Import the unittest classes.
import TestVariantContext
import TestGvcfVariant
import TestVariantContextFile


# Create the the test loader and suite.
loader = unittest.TestLoader()
suite = unittest.TestSuite()


# Add the three test classes that test the functionality of each program class.
suite.addTests(loader.loadTestsFromModule(TestVariantContext))
suite.addTests(loader.loadTestsFromModule(TestGvcfVariant))
suite.addTests(loader.loadTestsFromModule(TestVariantContextFile))


# Run the tests.
runner = unittest.TextTestRunner(verbosity=3)
result = runner.run(suite)