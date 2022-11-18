#!/usr/bin/env python
import argparse
from cyvcf2 import VCF, Writer, Variant
from gtfparse import read_gtf
from pysam import FastaFile


class ClonalHaploidPhaser:

    def __init__(self, input_vcf, output_vcf, input_gtf, input_fasta):
        self.vcf_reader = VCF(input_vcf)
        self.vcf_writer = Writer(output_vcf, self.vcf_reader)
        self.gtf = read_gtf(input_gtf)
        self.cds_regions = self.gtf.loc[self.gtf.feature == "CDS", :]
        self.cds_regions.loc[:, "uid"] = self.cds_regions.loc[:, "start":"end"].apply(lambda x: "{}:{}".format(x[0], x[1]), axis=1)
        self.fasta = FastaFile(input_fasta)

    def _overlap_amino_acid(self, first_variant, first_overlapping_cds, second_variant, second_overlapping_cds) -> bool:
        overlap = False
        if first_variant is not None and second_variant is not None:
            shared_cds = set(first_overlapping_cds.keys()).intersection(set(second_overlapping_cds.keys()))
            for c in shared_cds:
                start, end = first_overlapping_cds.get(c)
                first_amino_acid_index = int((first_variant.POS - start) / 3)
                second_amino_acid_index = int((second_variant.POS - start) / 3)
                if first_amino_acid_index == second_amino_acid_index:
                    overlap = True
        return overlap

    def _get_middle_sequence(self, first_variant, second_variant):
        # NOTE: assumes single chromosome
        middle_sequence = ""
        offset = second_variant.POS - (first_variant.POS + len(first_variant.REF) - 1)
        if offset > 1:
            middle_sequence = self.fasta.fetch(
                first_variant.CHROM, first_variant.POS, first_variant.POS - 1 + offset)
        return middle_sequence

    def _merge_variants(self, first_variant, second_variant):
        # consider that variants may not be in consecutive positions... we need the reference genome for this
        merged_variant = first_variant
        middle_sequence = self._get_middle_sequence(first_variant, second_variant)
        merged_variant.REF = merged_variant.REF + middle_sequence + second_variant.REF
        merged_variant.ALT = [merged_variant.ALT[0] + middle_sequence + second_variant.ALT[0]]
        return first_variant

    def run(self):
        previous_overlapping_cds = {}
        previous_variant = None
        variants_buffer = []

        variant: Variant
        for variant in self.vcf_reader:
            if variant.FILTER is None and variant.INFO.get('vafator_af', 1.0) >= 0.8:
                position = variant.POS
                overlapping_cds = {cds.uid: (cds.start, cds.end) for _, cds in
                                   self.cds_regions[
                                       (self.cds_regions.start <= position) &
                                       (self.cds_regions.end >= position)].iterrows()}
                if self._overlap_amino_acid(variant, overlapping_cds, previous_variant, previous_overlapping_cds):
                    # merge variants
                    previous_variant = self._merge_variants(first_variant=previous_variant, second_variant=variant)
                else:
                    # write previous variant
                    if previous_variant is not None:
                        variants_buffer.append(previous_variant)
                    previous_variant = variant
                previous_overlapping_cds = overlapping_cds
            else:
                # subclonal variants are not merged, this would require looking into read support
                variants_buffer.append(variant)

        # writes last variant
        if previous_variant:
            variants_buffer.append(previous_variant)

        # sorts the variants by position before writing them all
        for variant in sorted(variants_buffer, key=lambda v: v.POS):
            try:
                self.vcf_writer.write_record(variant)
            except Exception as e:
                pass

        self.vcf_reader.close()
        self.vcf_writer.close()


def main():
    parser = argparse.ArgumentParser(description="Joins PASS variants happening in same amino acid.")

    parser.add_argument("-i", "--input-vcf",
                        dest="vcf_in",
                        required=True,
                        help="Input vcf file listing somatic variants (expects sorted VCF).")
    parser.add_argument("-g", "--gtf",
                        dest="gtf",
                        required=True,
                        help="Input GTF defining the coding regions.")
    parser.add_argument("-f", "--fasta",
                        dest="fasta",
                        required=True,
                        help="Input FASTA reference genome.")
    parser.add_argument("-o", "--output-vcf",
                        dest="vcf_out",
                        required=True,
                        help="Output vcf file")

    args = parser.parse_args()

    ClonalHaploidPhaser(
        input_vcf=args.vcf_in, output_vcf=args.vcf_out, input_gtf=args.gtf, input_fasta=args.fasta).run()


if __name__ == "__main__":
    main()
