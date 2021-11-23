from unittest import TestCase
from .reference_genome import ReferenceGenomeReader


class TestReferenceGenomeReader(TestCase):

    def test_reference_genome_reader(self):
        reader = ReferenceGenomeReader(
            fasta_file_path="./reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
        base = reader.read_position_from_reference(chromosome="MN908947.3", position=123)
        self.assertEqual(base, "C")
        base = reader.read_position_from_reference(chromosome="MN908947.3", position=123, length=2)
        self.assertEqual(base, "CG")
        base = reader.read_position_from_reference(chromosome="MN908947.3", position=123, length=5)
        self.assertEqual(base, "CGCAG")
