from unittest import TestCase
from .assembly_variant_caller import AssemblyVariantCaller


class TestCountryParser(TestCase):

    def test_assembly_variant_caller(self):
        caller = AssemblyVariantCaller(
            match_score=2, mismatch_score=-1, open_gap_score=-3, extend_gap_score=-0.1)
        # no mutations
        variants = caller.call_variants(sequence="ACGTACGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)
        # SNV
        variants = caller.call_variants(sequence="ACGTCCGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "A")
        self.assertEqual(snv.alternate, "C")
        self.assertEqual(snv.position, 4)
        # deletion
        variants = caller.call_variants(
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG",
            sequence="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGCAAAGTGACG")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "TGGA")
        self.assertEqual(snv.alternate, "T")
        self.assertEqual(snv.position, 44)
        # insertion
        variants = caller.call_variants(
            sequence= "CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGCCCGAGGCAAAGTGACG",
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG")
        self.assertEqual(len(variants), 1)
        snv = variants[0]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GCCC")
        self.assertEqual(snv.position, 45)
        # another insertion
        variants = caller.call_variants(
            sequence= "CTGGTGTGAGTCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGCCCGAGGCAAAGTGACG",
            reference="CTGGTGTGAGCCTGGTCACCAGGGTGGTAGGACAGACCCTCCTCTGGAGGCAAAGTGACG")
        self.assertEqual(len(variants), 2)
        snv = variants[1]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GCCC")
        self.assertEqual(snv.position, 45)
        snv = variants[0]
        self.assertEqual(snv.reference, "G")
        self.assertEqual(snv.alternate, "GT")
        self.assertEqual(snv.position, 9)

    def test_ambiguous_bases(self):
        caller = AssemblyVariantCaller(
            match_score=2, mismatch_score=-1, open_gap_score=-3, extend_gap_score=-0.1)
        # no mutations
        variants = caller.call_variants(sequence="ACGTACGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)
        # ambiguous calls
        variants = caller.call_variants(sequence="ACGTWCGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)
        variants = caller.call_variants(sequence="ACGTNCGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)
        variants = caller.call_variants(sequence="ACGTAAAZRCCCCGT", reference="ACGTACGT")
        self.assertEqual(len(variants), 0)


