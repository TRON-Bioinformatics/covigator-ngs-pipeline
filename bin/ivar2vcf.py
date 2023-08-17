#!/usr/bin/env python3
import os
from argparse import ArgumentParser
from typing import List
import pandas as pd
from dataclasses import dataclass
from reference_genome import ReferenceGenomeReader


@dataclass
class Variant:
    chromosome: str
    position: int
    reference: str
    alternate: str
    passed: bool

    def to_vcf_line(self):
        # transform 0-based position to 1-based position
        return [self.chromosome, str(self.position), ".", self.reference, self.alternate, ".", "PASS" if self.passed else "LOW_QUALITY", "."]


class Ivar2Vcf:

    def __init__(self, ivar, reference_genome_reader: ReferenceGenomeReader):
        self.ivar_df = pd.read_csv(ivar, sep="\t")
        self.reference_genome_reader = reference_genome_reader
        self.chromosome = self.reference_genome_reader.chromosome

    def _is_ambiguous_base(self, base):
        return base not in "ACGT"

    def parse(self) -> List[Variant]:
        variants = []
        # remove repeated mutations
        self.ivar_df = self.ivar_df[["REGION", "POS", "REF", "ALT", "PASS"]].drop_duplicates()
        for index, row in self.ivar_df.iterrows():
            if len(row["ALT"]) == 1:
                # SNV
                variants.append(
                    Variant(
                        chromosome=self.chromosome,
                        position=row["POS"],
                        reference=row["REF"],
                        alternate=row["ALT"],
                        passed=bool(row["PASS"])))
            elif "+" in row["ALT"]:
                # insertion
                variants.append(
                    Variant(
                        chromosome=self.chromosome,
                        position=row["POS"],
                        reference=row["REF"],
                        alternate=row["REF"] + row["ALT"][1:],
                        passed=bool(row["PASS"])))
            elif "-" in row["ALT"]:
                # deletion
                deleted_sequence = self.reference_genome_reader.read_position_from_reference(
                    position=row["POS"] + 1, length=len(row["ALT"]) - 1)
                variants.append(
                    Variant(
                        chromosome=self.chromosome,
                        position=row["POS"],
                        reference=row["REF"] + deleted_sequence,
                        alternate=row["REF"],
                        passed=bool(row["PASS"])))
        return variants


def write_vcf(mutations, output_vcf, chromosome):
    with open(output_vcf, "w") as vcf_out:
        header = (
            "##fileformat=VCFv4.0",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">",
            "##FILTER=<ID=LOW_QUALITY,Description=\"Low quality variant call\">",
            "##contig=<ID={chromosome}>".format(chromosome=chromosome),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        )
        for row in header:
            vcf_out.write(row + "\n")
        for mutation in mutations:
            vcf_out.write("\t".join(mutation.to_vcf_line()) + "\n")


def main():
    parser = ArgumentParser(description="Ivar2VCF")
    parser.add_argument("--ivar", dest="ivar", help="The Ivar file to be converted into a VCF", required=True)
    parser.add_argument("--fasta", dest="fasta", help="The FASTA reference genome", required=True)
    parser.add_argument("--output-vcf", dest="output_vcf", help="The output VCF", required=True)

    args = parser.parse_args()

    assert os.path.exists(args.ivar), "Ivar file {} does not exist!".format(args.ivar)

    reference_genome_reader = ReferenceGenomeReader(fasta_file_path=args.fasta)
    parser = Ivar2Vcf(ivar=args.ivar, reference_genome_reader=reference_genome_reader)
    variants = parser.parse()
    write_vcf(mutations=variants, output_vcf=args.output_vcf, chromosome=reference_genome_reader.chromosome)


if __name__ == '__main__':
    main()
