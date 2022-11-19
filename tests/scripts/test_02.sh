#!/bin/bash

##################################################################################
# FASTQ input
# single-end reads
# --keep_intermediate has BAM files in output
##################################################################################
echo "Running CoVigator pipeline test 2"
source bin/assert.sh
output=tests/output/test2
nextflow main.nf -profile test,conda --name test_data \
--output $output \
--fastq1 tests/test_data/test_data_1.fastq.gz --keep_intermediate \
--skip_ivar --skip_bcftools --skip_gatk

test -s $output/test_data.lofreq.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.json || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.html || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.coverage.tsv || { echo "Missing coverage output file!"; exit 1; }
test -s $output/test_data.depth.tsv || { echo "Missing depth output file!"; exit 1; }
test -s $output/test_data.lofreq.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }
test -s $output/test_data.lofreq.fasta || { echo "Missing FASTA output file!"; exit 1; }

# these are the intermediate files kept by --keep_intermediate
test -s $output/test_data.preprocessed.bam || { echo "Missing BAM file!"; exit 1; }
test -s $output/test_data.preprocessed.bai || { echo "Missing BAI file!"; exit 1; }

assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | wc -l` 11 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep PASS | wc -l` 2 "Wrong number of variants"

assert_eq `cat $output/test_data.lofreq.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
