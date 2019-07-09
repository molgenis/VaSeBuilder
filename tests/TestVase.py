import unittest
import logging

from vase import VaSe
from ParamChecker import ParamChecker


class TestVase(unittest.TestCase):
    def setUp(self):
        self.vase = VaSe()
        self.paramchecker = ParamChecker()
        self.log_out_location = "testdata/outDir"
        self.valid_varlist_file_loc = "testdata/variantlist.txt"
        self.valid_config_file_loc = "testdata/configdata.cfg"
        self.invalid_varlist_file_loc = "testdata/invalid_variantlist.txt"
        self.invalid_config_file_loc = "tesdata/invalid_configdata.cfg"

    # Tests that the logger is started
    def test_start_logger(self):
        vase_logger = self.vase.start_logger(self.paramchecker, self.log_out_location)
        self.assertIsInstance(vase_logger, logging.Logger, "A logging.Logger object should have been returned")

    # Tests that the logging level of the logger is indeed set to DEBUG
    def test_start_logger_isindebug(self):
        is_debug = True
        vase_logger = self.vase.start_logger(self.paramchecker, self.log_out_location, is_debug)
        self.assertTrue(vase_logger.getEffectiveLevel() == 10, "The log level should have been set to DEBUG")

    # Tests that the variant list is read correctly
    def test_read_variant_list(self):
        varlist_answer = {}
        self.assertDictEqual(self.vase.read_variant_list(self.valid_varlist_file_loc), varlist_answer,
                             f"The read variant list should have been {varlist_answer}")

    # Tests that an empty variant list is returned if no file is submitted
    def test_read_variant_list_nofile(self):
        varlist_answer = {}
        non_existing_varlist = ""
        self.assertDictEqual(self.vase.read_variant_list(non_existing_varlist), varlist_answer,
                             "Both variant lists should have been empty")

    # Tests that an exception is raised if an invalid variant list file is provided
    def test_read_variant_list_invalidfile(self):
        varlist_answer = {}
        self.assertDictEqual(self.vase.read_config_file(self.invalid_varlist_file_loc), varlist_answer,
                             "The read variant list should have been empty")

    # Tests that a valid config file is read correctly.
    def test_read_config_file(self):
        config_data_answer = {"donorvcf": "testdata/vcfDir/vcflistfile.txt",
                              "donorbam": "testdata/bamDir/bamlistfile.txt",
                              "acceptorbam": "testdata/valbam/SRR1039513.bam",
                              "templatefq1": "testdata/fqDir/SRR1039513_1.fq.gz",
                              "templatefq2": "testdata/fqDir/SRR1039513_2.fq.gz",
                              "out": "testdata/outDir",
                              "reference": "testdata/ref/reference.fa"}
        self.assertDictEqual(self.vase.read_config_file(self.valid_config_file_loc), config_data_answer,
                             f"The read config data should have been: {config_data_answer}")

    # Tests that an empty config data map is returned when no file is supplied
    def test_read_config_file_nofile(self):
        config_data_answer = {}
        non_existing_config_file = ""
        self.assertDictEqual(self.vase.read_config_file(non_existing_config_file), config_data_answer,
                             "Both config data maps should have been empty")

    # Tests that an empty config data map is returned when a config file with an incorrect format is supplied
    def test_read_config_file_invalidfile(self):
        config_data_answer = {}
        self.assertDictEqual(self.vase.read_config_file(self.invalid_config_file_loc), config_data_answer,
                             "Both config data maps should have been empty")
