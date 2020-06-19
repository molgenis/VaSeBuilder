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
        self.valid_filter_header = tuple(["Aap", "Noot", "Mies", "Wim", "Does"])
        self.valid_header_line = "Aap\tNoot\tMies\tWim\tDoes\n"

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

    # Tests that a donorfastq list file is read correctly
    def test_read_donor_fastq_list_file(self):
        donorfqlist_answer = [["aR1.fq", "aR2.fq"], ["bR1.fq", "bR2.fq"], ["cR1.fq", "cR2.fq"]]
        self.assertListEqual(self.vase.read_donor_fastq_list_file("donorfqlist.txt"), donorfqlist_answer,
                             f"The rad donor fastq list file should have been {donorfqlist_answer}")

    # Tests that an empty donor fastq list is returned when no file is supplied
    def test_read_donor_fastq_list_file_nofile(self):
        donorfqlist_answer = []
        self.assertListEqual(self.vase.read_donor_fastq_list_file("nofile.txt"), donorfqlist_answer,
                             "The donor fastq list should have been empty")

    # Tests that a header is indeed valid.
    def test_is_valid_header(self):
        header_line = "Aap\tNoot\tMies\tWim\tDoes\n"
        self.assertTrue(self.vase.is_valid_header(header_line, self.valid_filter_header),
                        "The header line should have been correct")

    # Tests that an all CAPS header line with the correct fields is still found to be correct.
    def test_is_valid_header_capsline(self):
        header_line = "AAP\tNOOT\tMIES\tWIM\tDOES\n"
        self.assertTrue(self.vase.is_valid_header(header_line, self.valid_filter_header),
                        f"The all caps header line {header_line} should have been correct.")

    # Tests that a complete lowercase header line with the correct fields is still found to be correct.
    def test_is_valid_header_lowline(self):
        header_line = "aap\tnoot\tmies\twim\tdoes\n"
        self.assertTrue(self.vase.is_valid_header(header_line, self.valid_filter_header),
                        f"The lower case header line {header_line} should have been correct.")

    # Tests that a header line with incorrect fields is indeed found to be incorrect.
    def test_is_valid_header_invalidline(self):
        header_line = "does\twim\mies\tnoot\taap\n"
        self.assertFalse(self.vase.is_valid_header(header_line, self.valid_filter_header),
                         f"The header line {header_line} should have been incorrect.")

    # Tests that a filter that should be in the header is indeed found to be in the header.
    def test_filter_in_header(self):
        filter_name = "Does"
        self.assertTrue(self.vase.filter_in_header(filter_name, self.valid_header_line),
                        f"Filter name {filter_name} should have been in the header line {self.valid_header_line}")

    # Tests that a valid filter name in lowercase is still found to be in the header line.
    def test_filter_in_header_lowercase(self):
        filter_name = "does"
        self.assertTrue(self.vase.filter_in_header(filter_name, self.valid_header_line),
                        f"Lowercase filter name {filter_name} should have been in the header line"
                        f"{self.valid_header_line}")

    # Tests that a valid filter name in uppercase is still found to be in the header line.
    def test_filter_in_header_caps(self):
        filter_name = "DOES"
        self.assertTrue(self.vase.filter_in_header(filter_name, self.valid_header_line),
                        f"Uppercase filter name {filter_name} should have been in the header line"
                        f"{self.valid_header_line}")

    # Tests that an invalid filter name is indeed not in
    def test_filter_in_header_invalid(self):
        filter_name = "Piet"
        self.assertFalse(self.vase.filter_in_header(filter_name, self.valid_header_line),
                         f"Filter name {filter_name} should not have been in the header line "
                         f"{self.valid_header_line}.")

    # Tests that a valid filter name is at the expected index
    def test_get_filter_header_pos(self):
        filter_name = "Wim"
        filter_pos = 3
        self.assertEquals(self.vase.get_filter_header_pos(filter_name, self.valid_header_line), filter_pos,
                          f"The filter name {filter_name} should have been at index index {filter_pos}")

    # Tests that an invalid filter name is not in the header line and
    def test_get_filter_header_pos_none(self):
        filter_name = "Piet"
        self.assertIsNone(self.vase.get_filter_header_pos(filter_name, self.valid_header_line),
                          f"Filter name {filter_name} should not have been in the header line and therefore return a "
                          f"None.")

    # Tests that the correct priority index is returned.
    def test_determine_priority_index(self):
        priority_values = ("PATHOGENIC", "LIKELY PATHOGENIC", "BENIGN")
        filter_value = "PATHOGENIC"
        index_answer = 0
        self.assertEqual(self.vase.determine_priority_index(priority_values, filter_value), index_answer,
                         f"The returned priority level for {filter_value} should have been {index_answer}")

    # Tests that None is returned when a value is not within the filter.
    def test_determine_priority_index_none(self):
        priority_values = ("PATHOGENIC", "LIKELY PATHOGENIC", "BENIGN")
        filter_value = "AAP"
        self.assertIsNone(self.vase.determine_priority_index(priority_values, filter_value),
                          "The returned priority index should have been None")

    # Tests that None is returned when the filter is None.
    def test_determine_priority_index_nonefilter(self):
        priority_values = None
        filter_value = "PATHOGENIC"
        self.assertIsNone(self.vase.determine_priority_index(priority_values, filter_value),
                          "The returned priority index should have been None")

    # Tests that None is returned when a filter but no value is provided
    def test_determine_priority_index_nonevalue(self):
        priority_values = ("PATHOGENIC", "LIKELY PATHOGENIC", "BENIGN")
        filter_value = None
        self.assertIsNone(self.vase.determine_priority_index(priority_values, filter_value),
                          "The returned priority index should have been None")
