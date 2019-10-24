import unittest

from VaSeUtilRunner import VaSeUtilRunner


class TestVaSeUtilRunner(unittest.TestCase):
    def setUp(self):
        self.vsur = VaSeUtilRunner()

    def test_generate_config_file_from_command_shortnotation(self):
        short_command = "python vase.py -m F -v /testdata/donorvcfs.txt -b /testdata/donorbams.txt " \
                        "\-a /testdata/acceptorbam.bam " \
                        "-1 /testdata/template_L1_R1.fq.gz /testdata/template_L2_R1.fq.gz " \
                        "-2 /testdata/template_L1_R2.fq.gz /testdata/template_L2_R2.fq.gz -o /testdata/outdir " \
                        "-r /testdata/genomereference.fasta -l /testdata/outdir/out.log"

    def test_generate_config_file_from_command_longnotation(self):
        long_command = "python vase.py --runmode F --donorvcf /testdata/donorvcfs.txt" \
                       "--donorbam /testdata/donorbams.txt --acceptorbam /testdata/acceptorbam.bam" \
                       "--templatefq1 /testdata/template_L1_R1.fq.gz /testdata/template_L2_R1.fq.gz" \
                       "--templatefq2 /testdata/template_L1_R2.fq.gz /testdata/template_L2_R2.fq.gz" \
                       "--out /testdata/outdir --reference /testdata/genomereference.fasta" \
                       "--log /testdata/outdir/out.log"

    def test_generate_config_file_from_command_mixednotation(self):
        mixed_command = "python vase.py --runmode F -v /testdata/ -d /testdata/ -a /testdata/" \
                        "--templatefq1 /testdata/template_L1_R1.fq.gz /testdata/template_L2_R1.fq.gz" \
                        "--templatefq2 /testdata/template_L1_R2.fq.gz /testdata/template_L2_R2.fq.gz" \
                        "-o /testdata/outdir -r /testdata/genomereference.fasta --log /testdata/outdir/out.log"

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
