#!/usr/bin/env python
# Import required modules/libraries 
import unittest
import os
import sys

# Import required class
sys.path.insert(0, './../')
from ParamChecker import ParamChecker


# Unittest class for the ParamChecker test.
class TestParamChecker(unittest.TestCase):
    # Set up required 
    def setUp(self):
        self.param_check = ParamChecker()
        self.param_list = {"donorvcf": "testdata/vcfDir/vcflistfile.txt",
                           "donorbam": "testdata/bamDir/bamlistfile.txt",
                           "templatebam": "testdata/valbam/SRR1039513.bam",
                           "templatefq1": "testdata/fqDir/SRR1039513_1.fastq.gz",
                           "templatefq2": "testdata/fqDir/SRR1039513_2.fastq.gz",
                           "out": "testdata/outDir",
                           "reference": "testdata/ref/reference.fa",
                           'fastqout': "testdata/outDir",
                           'varcon': "testdata/outDir/varcon.txt",
                           "variantlist": "testdata/variantlist.txt"
                           }
        self.default_log_answer = "VaSeBuilder.log"
        self.default_varcon_answer = "varcon.txt"

    # Tests that the --log will be written to the log file.
    def test_check_log_pos(self):
        log_pos_answer = "testdata/outDir/vaseutest.log"
        self.assertEqual(self.param_check.check_log(log_pos_answer), log_pos_answer,
                         f"The log location should have been {log_pos_answer}")

    # Tests that the --log parameter has not been set and the log will therfore be written to VaSeBuilder.log 
    def test_check_log_no_param_set(self):
        self.assertEqual(self.param_check.check_log(None), self.default_log_answer,
                         f"The log file location should have been {self.default_log_answer}")

    # Tests that the --log parameter has been set but the filename does not end with .log and will therefore be written
    # to VaSeBuilder.log
    def test_check_log_no_log(self):
        invalid_log_name = "testdata/outDir/logfile.file"
        self.assertEqual(self.param_check.check_log(invalid_log_name), self.default_log_answer,
                         f"The log location should have been {self.default_log_answer}")

    # Test that the --log parameter has been set but the filename does not end with .txt  and will therefore be written
    # to VaSeBuilder.log
    def test_check_log_no_txt(self):
        invalid_log_name = "testdata/outDir/logfile.file"
        self.assertEqual(self.param_check.check_log(invalid_log_name), self.default_log_answer,
                         f"The log location should have been {self.default_log_answer}")

    # Tests that the povided testdata folder 'bamDir' indeed contains two BAM files.
    def test_check_folder_contents_pos(self):
        self.assertEqual(self.param_check.check_folder_contents("testdata/bamDir", "bam"), 2,
                         "Two bam files should have been found")

    # Tests that the provided testdata folder 'bamDir' contains no VCF files
    def test_check_folder_contents_neg(self):
        self.assertEqual(self.param_check.check_folder_contents("testdata/bamDir", "vcf"), 0,
                         "No VCF files should have been found")

    # Tests that the provided testdata folder 'bamDir' exists and contains BAM files.
    def test_check_folders_exist_pos(self):
        self.assertListEqual(self.param_check.check_folders_exist(['testdata/bamDir'], 'bam'), ['testdata/bamDir'],
                             "The folder should have been testdata/bamDir")

    # Tests that the provided testdata folder 'bamDir' exists and does not contain VCF files.
    def test_check_folders_exist_neg(self):
        self.assertListEqual(self.param_check.check_folders_exist(['testdata/bamDir'], 'vcf'), [],
                             "There should be no folder")

    # Tests that a provided file indeed exists.
    def test_check_file_exists_pos(self):
        existing_file = "testdata/bamDir/SRR1039508.bam"
        self.assertTrue(self.param_check.check_file_exists(existing_file), f"{existing_file} should exist")

    # Tests whether a non existing file indeed does not exist.
    def test_check_file_exists_neg(self):
        nonexisting_file = ""
        self.assertFalse(self.param_check.check_file_exists(nonexisting_file),
                         f"Non existing file {nonexisting_file} should not exist")

    # Test that the output location for out.txt is valid.
    def test_is_valid_output_location_pos(self):
        valid_out_loc = "testdata/outDir/out.txt"
        self.assertTrue(self.param_check.is_valid_output_location(valid_out_loc),
                        f"{valid_out_loc} should have been a valid output location")

    # Test that the output location for out.txt is invalid.
    def test_is_valid_output_location_neg(self):
        invalid_out_loc = "testdata/doesnotexist/out.txt"
        self.assertFalse(self.param_check.is_valid_output_location(invalid_out_loc),
                         f"{invalid_out_loc} should have been an invalid output location")

    # Tests that the foldername a file is located in gets returned properly.
    def test_get_folder_name_pos(self):
        valid_path = "testdata/noSampleDir/noSampleBam.bam"
        valid_path_answer = "testdata/noSampleDir"
        self.assertEqual(self.param_check.get_folder_name(valid_path), valid_path_answer,
                         f"Folder name should have been {valid_path_answer}")
    
    # Test that the foldername is empty.
    def test_get_folder_name_neg(self):
        self.assertEqual(self.param_check.get_folder_name(""), "", "Folder names should both be empty.")

    # Tests that all provided parameters are ok.
    # First test should return True, all others False due to one parameter missing
    def test_check_parameters_pos(self):
        self.assertTrue(self.param_check.check_parameters(self.param_list))

    def test_check_parameters_no_donor_vcf(self):
        par_list = self.param_list.copy()
        par_list['donorvcf'] = 'testdata/doesnotexist.txt'    # Set donorvcf parameter tp nonexisting vcf list file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_donor_bam(self):
        par_list = self.param_list.copy()
        par_list['donorbam'] = 'testdata/doesnotexist.txt'    # Set donorbam parameter to nonexisting bam list file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_bam(self):
        par_list = self.param_list.copy()
        par_list["acceptorbam"] = "testdata/doesnotexist.bam"    # Set acceptorbam parameter to nonexisting bam file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_fq1(self):
        par_list = self.param_list.copy()
        par_list["templatefq1"] = "testdata/doesnotexist.fq"    # Set templatefq1 parameter to nonexisting fastq file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_acceptor_fq2(self):
        par_list = self.param_list.copy()
        par_list['templatefq2'] = 'testdata/doesnotexist.fq'    # Set templatefq2 parameter to nonexisting fastq file.
        self.assertFalse(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_fq_out(self):
        par_list = self.param_list.copy()
        par_list["fastqout"] = "testdata/doesnotexist/faqOut.fq"    # Set fastqout parameter to nonexisting location
        self.assertTrue(self.param_check.check_parameters(par_list))

    def test_check_parameters_no_varcon(self):
        par_list = self.param_list.copy()
        par_list["varcon"] = "testdata/doesnotexist/varcon.txt"    # Set varcon parameter to nonexisting location.
        self.assertTrue(self.param_check.check_parameters(par_list))

    def test_
