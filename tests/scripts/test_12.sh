#!/bin/bash

##################################################################################
# VCF input
##################################################################################
echo "Running CoVigator pipeline test 12"
source bin/assert.sh
output=tests/output/test12
nextflow main.nf -profile test,conda --name test_data \
	--output $output \
	--vcf tests/test_data/test_data.lofreq.vcf --keep_intermediate \
	--skip_sarscov2_annotations \
	--skip_pangolin \
	--skip_normalization

test -s $output/test_data.input.vcf.gz || { echo "Missing VCF output file!"; exit 1; }

assert_eq `zcat $output/test_data.input.vcf.gz | grep -v '#' | wc -l` 60 "Wrong number of variants"
# there is one mutation missing vafator_af as these annotations are missing in the input
assert_eq `zcat $output/test_data.input.vcf.gz | grep -v '#' | grep 'vafator_af' | wc -l` 59 "Wrong number of variants"
assert_eq `zcat $output/test_data.input.vcf.gz | grep -v '#' | grep PASS | wc -l` 6 "Wrong number of variants"
