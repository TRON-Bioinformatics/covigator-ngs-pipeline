#!/bin/bash

##################################################################################
# FASTA input
##################################################################################
echo "Running CoVigator pipeline test 4"
source bin/assert.sh
output=output/test4
nextflow main.nf -profile test,conda --name hCoV-19_NTXX \
	--output $output \
	--fasta test_data/hCoV-19_NTXX.fasta

test -s $output/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/hCoV-19_NTXX.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }

assert_eq `zcat $output/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 13 "Wrong number of variants"
assert_eq `cat $output/ERR4145453.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"