#!/usr/bin/env python
import argparse
import gzip


def get_header(vcf, sample):
    header = "\n".join(vcf._header_lines)
    if "##FORMAT=<ID=GT" not in header:
        header = header + '\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
    return header + "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sample)


def run(input_vcf, output_vcf, sample):
    vcf_lines = open(input_vcf, mode="r").readlines()

    header = [h for h in vcf_lines if h.startswith("##")]
    header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sample))

    variants = [v for v in vcf_lines if not v.startswith("#")]
    modified_variants = map(lambda v: v[0:-2] + "\tGT\t0/1\n", variants)

    vcf_writer = open(output_vcf, "w")
    vcf_writer.writelines(header)
    vcf_writer.writelines(modified_variants)

    vcf_writer.close()


def main():
    parser = argparse.ArgumentParser(description="Adds a single sample with an haploid genotype to the VCF file. "
                                                 "Any previous samples or genotype information is removed.")

    parser.add_argument("-i", "--input-vcf",
                        dest="vcf_in",
                        required=True,
                        help="Input vcf file listing somatic variants")
    parser.add_argument("-o", "--output-vcf",
                        dest="vcf_out",
                        required=True,
                        help="Output vcf file")
    parser.add_argument("-s", "--sample",
                        dest="sample",
                        default="sample",
                        help="Sample name default: sample).")

    args = parser.parse_args()

    run(args.vcf_in, args.vcf_out, args.sample)


if __name__ == "__main__":
    main()
