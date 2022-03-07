import os.path
from unittest import TestCase
from .phasing import ClonalHaploidPhaser
from cyvcf2 import VCF


class TestPhasing(TestCase):

    def test_reference_genome_reader(self):
        os.makedirs("./output/phasing", exist_ok=True)
        output_vcf = "./output/phasing/test_data.merged.vcf.gz"
        if os.path.exists(output_vcf):
            os.remove(output_vcf)
        ClonalHaploidPhaser(
            input_vcf="./test_data/test_data.lofreq.vcf",
            output_vcf=output_vcf,
            input_gtf="./reference/Sars_cov_2.ASM985889v3.101.gff3",
            input_fasta="./reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa"
        ).run()
        self.assertTrue(os.path.exists(output_vcf))
        self._assert_vcf(output_vcf, 56)

    def _assert_vcf(self, vcf, expected_count_variants):
        vcf_reader = VCF(vcf)
        count_variants = 0
        for v in vcf_reader:
            count_variants += 1
        self.assertEqual(count_variants, expected_count_variants)
