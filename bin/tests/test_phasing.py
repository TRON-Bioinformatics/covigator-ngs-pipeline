import os.path
from unittest import TestCase
from bin.phasing import ClonalHaploidPhaser
from cyvcf2 import VCF


class TestPhasing(TestCase):

    def test_phasing(self):
        os.makedirs(os.path.join(os.path.dirname(__file__), "output/phasing"), exist_ok=True)
        output_vcf = os.path.join(os.path.dirname(__file__), "output/phasing/test_data.merged.vcf.gz")
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        ClonalHaploidPhaser(
            input_vcf=os.path.join(os.path.dirname(__file__), "data/test_data.lofreq.vcf"),
            output_vcf=output_vcf,
            input_gtf=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.101.gff3"),
            input_fasta=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
        ).run()
        self.assertTrue(os.path.exists(output_vcf))
        self._assert_vcf(output_vcf, 60)

    def test_phasing_empty_vcf(self):
        os.makedirs(os.path.join(os.path.dirname(__file__), "output/phasing"), exist_ok=True)
        output_vcf = os.path.join(os.path.dirname(__file__), "output/phasing/test_data.merged.vcf.gz")
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        ClonalHaploidPhaser(
            input_vcf=os.path.join(os.path.dirname(__file__), "data/test_data.empty.vcf"),
            output_vcf=output_vcf,
            input_gtf=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.101.gff3"),
            input_fasta=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
        ).run()
        self.assertTrue(os.path.exists(output_vcf))
        self._assert_vcf(output_vcf, 0)

    def test_phasing_bug36(self):
        os.makedirs(os.path.join(os.path.dirname(__file__), "output/phasing"), exist_ok=True)
        output_vcf = os.path.join(os.path.dirname(__file__), "output/phasing/test.bug36.merged.vcf.gz")
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        ClonalHaploidPhaser(
            input_vcf=os.path.join(os.path.dirname(__file__), "data/test_data.bug36.vcf"),
            output_vcf=output_vcf,
            input_gtf=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.101.gff3"),
            input_fasta=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
        ).run()
        self.assertTrue(os.path.exists(output_vcf))
        self._assert_vcf(output_vcf, 2)

    def test_phasing_bug36_2_legitimate_overlapping_indels(self):
        os.makedirs(os.path.join(os.path.dirname(__file__), "output/phasing"), exist_ok=True)
        output_vcf = os.path.join(os.path.dirname(__file__), "output/phasing/test.bug36_2.merged.vcf.gz")
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        ClonalHaploidPhaser(
            input_vcf=os.path.join(os.path.dirname(__file__), "data/test_data.bug36_2.vcf"),
            output_vcf=output_vcf,
            input_gtf=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.101.gff3"),
            input_fasta=os.path.join(os.path.dirname(__file__), "../../reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa")
        ).run()
        self.assertTrue(os.path.exists(output_vcf))
        self._assert_vcf(output_vcf, 2)

    def _assert_vcf(self, vcf, expected_count_variants):
        vcf_reader = VCF(vcf)
        count_variants = 0
        for v in vcf_reader:
            count_variants += 1
        self.assertEqual(count_variants, expected_count_variants)
