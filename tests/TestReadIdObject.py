import unittest
from ReadIdObject import ReadIdObject


class TestReadIdObject(unittest.TestCase):
    def setUp(self):
        self.read_id_answer = "SEQ001"
        self.readidobj_answer = ReadIdObject(self.read_id_answer)

    def test_get_bam_read_id(self):
        self.assertEqual(self.readidobj_answer.get_bam_read_id(), self.read_id_answer, "Both read identifiers should "
                         f"have been {self.read_id_answer}")

    def test_set_bam_read_id(self):
        set_read_id_answer = "SEQ005"
        self.readidobj_answer.set_bam_read_id(set_read_id_answer)
        self.assertEqual(self.readidobj_answer.get_bam_read_id(), set_read_id_answer, "Both read identifiers should "
                         f"have been {set_read_id_answer}")
