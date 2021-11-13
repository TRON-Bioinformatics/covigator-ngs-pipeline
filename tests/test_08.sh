#!/bin/bash

##################################################################################
# FASTQ input
# paired-end reads
##################################################################################
echo "Running CoVigator pipeline test 8"
source bin/assert.sh
output=output/test8
nextflow main.nf -profile test,conda --name ERR4145453 \
	--output $output \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz --skip_sarscov2_annotations

test -s $output/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.ivar.tsv || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.json || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.fastp_stats.html || { echo "Missing VCF output file!"; exit 1; }
test -s $output/ERR4145453.coverage.tsv || { echo "Missing coverage output file!"; exit 1; }
test -s $output/ERR4145453.depth.tsv || { echo "Missing depth output file!"; exit 1; }
test -s $output/ERR4145453.depth.tsv || { echo "Missing deduplication metrics file!"; exit 1; }

assert_eq `zcat $output/ERR4145453.lofreq.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 2225 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.lofreq.normalized.annotated.vcf.gz | grep -v '#' | grep PASS | wc -l` 5 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.lofreq.normalized.annotated.vcf.gz | grep -v '#' | grep PFAM | wc -l` 0 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.bcftools.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 5 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.bcftools.normalized.annotated.vcf.gz | grep -v '#' | grep PFAM | wc -l` 0 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.gatk.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 6 "Wrong number of variants"
assert_eq `zcat $output/ERR4145453.gatk.normalized.annotated.vcf.gz | grep -v '#' | grep PFAM | wc -l` 0 "Wrong number of variants"