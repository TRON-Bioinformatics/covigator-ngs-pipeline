#!/bin/bash

output=output/test1
nextflow main.nf -profile test,conda --name ERR4145453 \
	--output $output \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz

test -s $output/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.json || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.html || { echo "Missing test 1 VCF output file!"; exit 1; }
test -s $output/ERR4145453.coverage.tsv || { echo "Missing test 1 coverage output file!"; exit 1; }
test -s $output/ERR4145453.depth.tsv || { echo "Missing test 1 depth output file!"; exit 1; }
