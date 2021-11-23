#!/usr/bin/env python
import os
from argparse import ArgumentParser
from typing import List
import pandas as pd
from dataclasses import dataclass
from reference_genome import ReferenceGenomeReader


DEFAULT_CHROMOSOME = "MN908947.3"


@dataclass
class Variant:
    position: int
    reference: str
    alternate: str
    passed: bool

    def to_vcf_line(self, chromosome):
        # transform 0-based position to 1-based position
        return [chromosome, str(self.position), ".", self.reference, self.alternate, ".", "PASS" if self.passed else "LOW_QUALITY", "."]


class Ivar2Vcf:

    def __init__(self, ivar, fasta, chromosome):
        self.ivar_df = pd.read_csv(ivar, sep="\t")
        self.reference_genome_reader = ReferenceGenomeReader(fasta_file_path=fasta)
        self.chromosome = chromosome

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
                    Variant(position=row["POS"], reference=row["REF"], alternate=row["ALT"], passed=bool(row["PASS"])))
            elif "+" in row["ALT"]:
                # insertion
                variants.append(
                    Variant(position=row["POS"], reference=row["REF"], alternate=row["REF"] + row["ALT"][1:],
                            passed=bool(row["PASS"])))
            elif "-" in row["ALT"]:
                # deletion
                deleted_sequence = self.reference_genome_reader.read_position_from_reference(
                    chromosome=self.chromosome, position=row["POS"] + 1, length=len(row["ALT"]) - 1)
                variants.append(
                    Variant(position=row["POS"], reference=row["REF"] + deleted_sequence, alternate=row["REF"],
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
            vcf_out.write("\t".join(mutation.to_vcf_line(chromosome=chromosome)) + "\n")


def main():
    parser = ArgumentParser(description="Ivar2VCF")
    parser.add_argument("--ivar", dest="ivar", help="The Ivar file to be converted into a VCF", required=True)
    parser.add_argument("--fasta", dest="fasta", help="The FASTA reference genome", required=True)
    parser.add_argument("--output-vcf", dest="output_vcf", help="The output VCF", required=True)
    parser.add_argument("--chromosome", dest="chromosome",
                        help="The chromosome to be used in the output VCF. Beware only one chromosome is supported!",
                        default=DEFAULT_CHROMOSOME)

    args = parser.parse_args()

    assert os.path.exists(args.ivar), "Ivar file {} does not exist!".format(args.ivar)

    parser = Ivar2Vcf(ivar=args.ivar, chromosome=args.chromosome, fasta=args.fasta)
    variants = parser.parse()
    write_vcf(mutations=variants, output_vcf=args.output_vcf, chromosome=args.chromosome)


if __name__ == '__main__':
    main()
