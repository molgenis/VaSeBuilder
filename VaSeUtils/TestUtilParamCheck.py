import unittest

from UtilParamCheck import UtilParamCheck


class TestUtilParamCheck(unittest.TestCase):
    def setUp(self):
        self.vupc = UtilParamCheck()

    def test_required_params_set_positive(self):
        utiltorun = "acceptorcheck"
        utilparams = {"util": "acceptorcheck", "varcon": "/testdata/variantcontext.txt",
                      "vasefq1": "/testdata/template_r1.fq", "vasefq2": "/testdata/template_r2.fq"}
        self.assertTrue(self.vupc.required_params_set(utiltorun, utilparams),
                        f"VaSeUtil {utiltorun} should have had all parameters set correctly.")

    def test_required_params_set_negative(self):
        utiltorun = "acceptorcheck"
        utilparams = {"util": "acceptorcheck", "varcon": "/testdata/variantcontext.txt",
                      "vasefq1": "/testdata/template_r1.fq"}
        self.assertFalse(self.vupc.required_params_set(utiltorun, utilparams),
                         f"VaSeUtil {utiltorun} should not have all parameters set correctly.")

    def test_get_required_parameters_positive(self):
        utiltorun = "donorcheck"
        required_params_answer = ["varcon", "vasefq1", "vasefq2"]
        self.assertListEqual(self.vupc.get_required_parameters(utiltorun), required_params_answer,
                             f"VaSeUtil {utiltorun} should have returned required parameters {required_params_answer}")

    def test_get_required_parameters_negative(self):
        utiltorun = "nonutil"
        self.assertIsNone(self.vupc.get_required_parameters(utiltorun),
                          f"The required util {utiltorun} should not exists and therefore return None")

    def test_get_not_set_parameters_positive(self):
    def test_get_not_set_parameters_negative(self):

    def test_get_set_optional_parameters_positive(self):
    def test_get_set_optional_parameters_negative(self):

    def test_get_unused_optional_parameters_positive(self):
    def test_get_unused_optional_parameters_negative(self):

    def test_get_valid_vaseutils(self):
        valid_vaseutils = ['acceptorcheck', 'acceptorreadinfo', 'checkdonorfiles', 'checkfastq', 'compareacceptor',
                           'comparedonor', 'comparefastq', 'comparevarcon', 'donorcheck', 'donorreadinfo', 'loginfo',
                           'subsetvarcon', 'subvcfvarcon', 'subvcfvarlist', 'unmappedinfo', 'varcondata' ,'vcfinvarcon']
        self.assertListEqual(self.vupc.get_valid_vaseutils(), valid_vaseutils,
                             f"The list valid VaSeUtils should have been {valid_vaseutils}")
