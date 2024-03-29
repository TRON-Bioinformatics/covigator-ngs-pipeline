#!/bin/bash

##################################################################################
# FASTQ input
# paired-end reads
##################################################################################
echo "Running CoVigator pipeline test 1"
source bin/assert.sh
output=tests/output/test1
nextflow main.nf -profile test,conda --name test_data \
	--output $output \
	--fastq1 tests/test_data/test_data_1.fastq.gz \
	--fastq2 tests/test_data/test_data_2.fastq.gz --keep_intermediate

test -s $output/test_data.bcftools.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.gatk.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.lofreq.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.lofreq.vcf.gz.tbi || { echo "Missing VCF index output file!"; exit 1; }
test -s $output/test_data.ivar.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.ivar.tsv || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.json || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.html || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.coverage.tsv || { echo "Missing coverage output file!"; exit 1; }
test -s $output/test_data.depth.tsv || { echo "Missing depth output file!"; exit 1; }
test -s $output/test_data.lofreq.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }
test -s $output/test_data.lofreq.fasta || { echo "Missing FASTA output file!"; exit 1; }
test -s $output/test_data.gatk.fasta || { echo "Missing FASTA output file!"; exit 1; }
test -s $output/test_data.bcftools.fasta || { echo "Missing FASTA output file!"; exit 1; }
test -s $output/test_data.ivar.fasta || { echo "Missing FASTA output file!"; exit 1; }

assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | wc -l` 54 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep 'vafator_af' | wc -l` 54 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep PASS | wc -l` 2 "Wrong number of variants"
assert_eq `zcat $output/test_data.bcftools.vcf.gz | grep -v '#' | wc -l` 13 "Wrong number of variants"
assert_eq `zcat $output/test_data.gatk.vcf.gz | grep -v '#' | wc -l` 11 "Wrong number of variants"

assert_eq `cat $output/test_data.lofreq.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
