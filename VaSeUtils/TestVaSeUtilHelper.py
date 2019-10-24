import unittest

from VaSeUtilHelper import VaSeUtilHelper


class TestVaSeUtilHelper(unittest.TestCase):
    def setUp(self):
        self.vsuh = VaSeUtilHelper()

    def test_passes_filter_positive(self):
        val_to_use = "aap"
        filter_to_use = ["aap", "noot", "mies"]
        self.assertTrue(self.vsuh.passes_filter(),
                        f"Value {val_to_use} should have passed filter {filter_to_use}")

    def test_passes_filter_negative(self):
        val_to_use = "does"
        filter_to_use = ["aap", "noot", "mies"]
        self.assertFalse(self.vsuh.passes_filter,
                         f"Value {val_to_use} should not have passed filter {filter_to_use}")

    def test_get_accdon_context_fields(self):
        fields_answer = {1: "Context ID",
                         2: "Sample ID",
                         3: "Context chrom",
                         4: "Context origin",
                         5: "Context start",
                         6: "Context end",
                         7: "Context reads"}
        self.assertDictEqual(self.vsuh.get_accdon_context_fields(), fields_answer,
                             f"The returned map of acceptor/donor context fields should have been {fields_answer}")

    def test_get_config_parameter_name_posshort(self):
        short_param_flag = "-a"
        parameter_name_answer = "ACCEPTORBAM"
        self.assertEqual(self.vsuh.get_config_parameter_name(short_param_flag), parameter_name_answer,
                         f"Parameter name for flag {short_param_flag} should have been {parameter_name_answer}")

    def test_get_config_parameter_name_poslong(self):
        long_param_flag = "--acceptorbam"
        parameter_name_answer = "ACCEPTORBAM"
        self.assertEqual(self.vsuh.get_config_parameter_name(long_param_flag), parameter_name_answer,
                         f"Parameter name for flag {long_param_flag} should have been {parameter_name_answer}")

    def test_get_config_parameter_name_negshort(self):
        short_param_flag = "-q"
        self.assertIsNone(self.vsuh.get_config_parameter_name(short_param_flag),
                          f"Parameter flag {short_param_flag} should not exist and therefore return None")

    def test_get_config_parameter_name_neglong(self):
        long_param_flag = "--quit"
        self.assertIsNone(self.vsuh.get_config_parameter_name(),
                          f"Parameter flag {long_param_flag} should not exist and therefore return None")

    def test_is_valid_parameter_flag_posshort(self):
        short_param_flag = "-m"
        self.assertTrue(self.vsuh.is_valid_parameter_flag(short_param_flag),
                        f"Parameter flag {short_param_flag} should have been a valid parameter flag.")

    def test_is_valid_parameter_flag_poslong(self):
        long_param_flag = "--runmode"
        self.assertTrue(self.vsuh.is_valid_parameter_flag(long_param_flag),
                        f"Parameter flag {long_param_flag} should have been a valid parameter flag.")

    def test_is_valid_parameter_flag_negshort(self):
        short_param_flag = "-q"
        self.assertFalse(self.vsuh.is_valid_parameter_flag(short_param_flag),
                         f"Parameter flag {short_param_flag} should have been an invalid parameter flag.")

    def test_is_valid_parameter_flag_neglong(self):
        long_param_flag = "--quit"
        self.assertFalse(self.vsuh.is_valid_parameter_flag(long_param_flag),
                         f"Parameter flag {long_param_flag} should have been an invalid parameter flag.")
