from unittest import TestCase
from bin.reference_genome import ReferenceGenomeReader
import os


class TestReferenceGenomeReader(TestCase):

    def test_reference_genome_reader(self):
        reader = ReferenceGenomeReader(
            fasta_file_path=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa"))
        base = reader.read_position_from_reference(position=123)
        self.assertEqual(base, "C")
        base = reader.read_position_from_reference(position=123, length=2)
        self.assertEqual(base, "CG")
        base = reader.read_position_from_reference(position=123, length=5)
        self.assertEqual(base, "CGCAG")
