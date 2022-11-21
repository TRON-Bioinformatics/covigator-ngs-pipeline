#!/bin/bash

##################################################################################
# FASTQ input
# paired-end reads
##################################################################################
echo "Running CoVigator pipeline test 1"
source bin/assert.sh
output=tests/output/test12
nextflow main.nf -profile test,conda --name test_data \
	--output $output \
	--vcf tests/test_data/test_data.lofreq.vcf --keep_intermediate \
	--skip_sarscov2_annotations \
	--skip_pangolin \
	--skip_normalization

test -s $output/test_data.input.vcf.gz || { echo "Missing VCF output file!"; exit 1; }

assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | wc -l` 54 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep 'vafator_af' | wc -l` 54 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep PASS | wc -l` 2 "Wrong number of variants"
