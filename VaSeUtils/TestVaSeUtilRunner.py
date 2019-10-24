import unittest

from VaSeUtilRunner import VaSeUtilRunner


class TestVaSeUtilRunner(unittest.TestCase):
    def setUp(self):
        self.vsur = VaSeUtilRunner()

    def test_is_multivalue_parameter_positive(self):
        multival_param = "templatefq1"
        self.assertTrue(self.vsur.is_multivalue_parameter(multival_param),
                        f"The parameter name {multival_param} should have been a multi-value parameter")

    def test_is_multivalue_parameter_negative(self):
        nonmultival_param = "debug"
        self.assertFalse(self.vsur.is_multivalue_parameter(nonmultival_param),
                         f"The parameter name {nonmultival_param} should not have been a multi-value parameter.")

    def test_is_nonvalue_parameter_positive(self):
        nonvalue_param = "debug"
        self.assertTrue(self.vsur.is_nonvalue_parameter(nonvalue_param),
                        f"The parameter name {nonvalue_param} should have been a non-value parameter")

    def test_is_nonvalue_parameter_negative(self):
        value_param = "acceptorbam"
        self.assertFalse(self.vsur.is_nonvalue_parameter(value_param),
                         f"The parameter name {value_param} should not have been a non-value parameter")
