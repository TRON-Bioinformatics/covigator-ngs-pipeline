#!/bin/bash

##################################################################################
# FASTA input
##################################################################################
echo "Running CoVigator pipeline test 4"
source bin/assert.sh
output=output/test4
nextflow main.nf -profile test,conda --name test_data \
	--output $output \
	--fasta test_data/test_data.fasta

test -s $output/test_data.assembly.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }

assert_eq `zcat $output/test_data.assembly.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 13 "Wrong number of variants"
assert_eq `cat $output/test_data.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"