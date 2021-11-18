#!/bin/bash

##################################################################################
# FASTQ input
# paired-end reads
# use --input_fastqs_list
##################################################################################
echo "Running CoVigator pipeline test 6"
source bin/assert.sh
output=output/test6
echo -e "ERR4145453\t"`pwd`"/test_data/ERR4145453_1.fastq.gz\t"`pwd`"/test_data/ERR4145453_2.fastq.gz\n" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --input_fastqs_list test_data/test_input.txt \
	--library paired --output $output

test -s $output/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.ivar.tsv || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.json || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.html || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.coverage.tsv || { echo "Missing coverage output file!"; exit 1; }
test -s $output/ERR4145453.depth.tsv || { echo "Missing depth output file!"; exit 1; }
test -s $output/ERR4145453.depth.tsv || { echo "Missing deduplication metrics file!"; exit 1; }
test -s $output/ERR4145453.bcftools.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }
test -s $output/ERR4145453.gatk.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }
test -s $output/ERR4145453.lofreq.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }

assert_eq `zcat $output/ERR4145453.lofreq.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 2225 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.lofreq.normalized.annotated.vcf.gz | grep -v '#' | grep PASS | wc -l` 5 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.bcftools.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 5 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.gatk.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 6 "Wrong number of variants"

assert_eq `cat $output/ERR4145453.gatk.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
assert_eq `cat $output/ERR4145453.bcftools.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
assert_eq `cat $output/ERR4145453.lofreq.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
