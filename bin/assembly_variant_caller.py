#!/usr/bin/env python
import os
from argparse import ArgumentParser
from dataclasses import dataclass
from Bio import Align, SeqIO
from Bio.Align import PairwiseAlignment
from typing import List

DEFAULT_EXTEND_GAP_SCORE = -0.1
DEFAULT_OPEN_GAP_SCORE = -3
DEFAULT_MISMATCH_SCORE = 1
DEFAULT_MATCH_SCORE = 2


@dataclass
class Variant:
    chromosome: str
    position: int
    reference: str
    alternate: str

    def to_vcf_line(self):
        # transform 0-based position to 1-based position
        return [self.chromosome, str(self.position + 1), ".", self.reference, self.alternate, ".", "PASS", "."]


class AssemblyVariantCaller:

    def __init__(self, match_score, mismatch_score, open_gap_score, extend_gap_score):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.open_gap_score = open_gap_score
        self.extend_gap_score = extend_gap_score

    def call_variants(self, sequence: str, reference: str, chromosome: str) -> List[Variant]:
        alignment = self._run_alignment(sequence=sequence, reference=reference)
        variants = self._call_mutations(alignment, chromosome=chromosome)
        return variants

    def _run_alignment(self, sequence: str, reference: str) -> PairwiseAlignment:
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match = self.match_score
        aligner.mismatch = self.mismatch_score
        aligner.open_gap_score = self.open_gap_score
        aligner.extend_gap_score = self.extend_gap_score
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        alignments = aligner.align(reference, sequence)
        return alignments[0]

    def _call_mutations(self, alignment: PairwiseAlignment, chromosome: str) -> List[Variant]:
        # CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        # MN908947.3      9924    .       C       T       228     .
        # DP=139;VDB=0.784386;SGB=-0.693147;RPB=0.696296;MQB=1;MQSB=1;BQB=0.740741;MQ0F=0;AC=1;AN=1;DP4=2,0,123,12;MQ=60
        # GT:PL   1:255,0
        alternate = alignment.query
        reference = alignment.target

        variants = []
        prev_ref_end = None
        prev_alt_end = None
        for (ref_start, ref_end), (alt_start, alt_end) in zip(alignment.aligned[0], alignment.aligned[1]):
            # calls indels
            # NOTE: it does not call indels at beginning and end of sequence
            if prev_ref_end is not None and prev_ref_end != ref_start:
                # deletion
                if ref_start - prev_ref_end <= 50:  # skips deletions longer than 50 bp
                    ref = str(reference[prev_ref_end - 1: ref_start])
                    if not any(self._is_ambiguous_base(r) for r in ref):  # do not call deletions with Ns
                        variants.append(Variant(
                            chromosome=chromosome,
                            position=prev_ref_end - 1,
                            reference=ref,
                            alternate=reference[prev_ref_end - 1]))
            elif prev_ref_end is not None and prev_alt_end != alt_start:
                # insertion
                if alt_start - prev_alt_end <= 50:  # skips insertions longer than 50 bp
                    ref = reference[prev_ref_end - 1]
                    alt = str(alternate[prev_alt_end:alt_start])
                    # do not call insertions with ambiguous bases
                    if not self._is_ambiguous_base(ref) and not any(self._is_ambiguous_base(a) for a in alt):
                        variants.append(Variant(
                            chromosome=chromosome,
                            position=prev_ref_end - 1,
                            reference=ref,
                            alternate=ref + alt))

            # calls SNVs
            for pos, ref, alt in zip(
                    range(ref_start, ref_end), reference[ref_start: ref_end], alternate[alt_start: alt_end]):
                # contiguous SNVs are reported separately
                # do not call insertions with ambiguous bases
                if ref != alt and not self._is_ambiguous_base(ref) and not self._is_ambiguous_base(alt):
                    variants.append(Variant(
                        chromosome=chromosome,
                        position=pos,
                        reference=ref,
                        alternate=alt))

            prev_ref_end = ref_end
            prev_alt_end = alt_end

        return variants

    def _is_ambiguous_base(self, base):
        return base not in "ACGT"


def write_vcf(mutations, output_vcf, chromosome):
    with open(output_vcf, "w") as vcf_out:
        header = (
            "##fileformat=VCFv4.0",
            "##FILTER=<ID=PASS,Description=\"All filters passed\">",
            "##contig=<ID={chromosome}>".format(chromosome=chromosome),
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        )
        for row in header:
            vcf_out.write(row + "\n")
        for mutation in mutations:
            vcf_out.write("\t".join(mutation.to_vcf_line()) + "\n")


def main():
    parser = ArgumentParser(description="Run Pipeline for testing")
    parser.add_argument("--fasta", dest="fasta",
                        help="The fasta file with the query sequence. Only one sequence is expected",
                        required=True)
    parser.add_argument("--reference", dest="reference",
                        help="The fasta file with the reference sequence. Only one sequence is expected",
                        required=True)
    parser.add_argument("--output-vcf", dest="output_vcf",
                        help="The path to the output VCF",
                        required=True)
    parser.add_argument("--match-score", dest="match_score", help="The score for a matching position",
                        default=DEFAULT_MATCH_SCORE)
    parser.add_argument("--mismatch-score", dest="mismatch_score", help="The score for a mismatching position",
                        default=-DEFAULT_MISMATCH_SCORE)
    parser.add_argument("--open-gap-score", dest="open_gap_score", help="The score for opening a gap",
                        default=DEFAULT_OPEN_GAP_SCORE)
    parser.add_argument("--extend-gap-score", dest="extend_gap_score", help="The score for extending a gap",
                        default=DEFAULT_EXTEND_GAP_SCORE)

    args = parser.parse_args()

    assert os.path.exists(args.fasta), "Fasta file {} does not exist!".format(args.fasta)
    assert os.path.exists(args.reference), "Fasta file {} does not exist!".format(args.reference)

    query = next(SeqIO.parse(args.fasta, "fasta"))
    reference = next(SeqIO.parse(args.reference, "fasta"))
    variant_caller = AssemblyVariantCaller(
        match_score=float(args.match_score),
        mismatch_score=float(args.mismatch_score),
        open_gap_score=float(args.open_gap_score),
        extend_gap_score=float(args.extend_gap_score)
    )
    variants = variant_caller.call_variants(sequence=query.seq, reference=reference.seq, chromosome=reference.id)
    write_vcf(mutations=variants, output_vcf=args.output_vcf, chromosome=reference.id)


if __name__ == '__main__':
    main()
