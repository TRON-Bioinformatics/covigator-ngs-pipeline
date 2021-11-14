#!/bin/bash

##################################################################################
# FASTA input
# use --input_fastas_list
##################################################################################
echo "Running CoVigator pipeline test 5"
source bin/assert.sh
output=output/test5
echo -e "hCoV-19_NTXX\t"`pwd`"/test_data/hCoV-19_NTXX.fasta\n" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --input_fastas_list test_data/test_input.txt \
	--output $output

test -s $output/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz || { echo "Missing VCF output file!"; exit 1; }

assert_eq `zcat $output/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz | grep -v '#' | wc -l` 13 "Wrong number of variants"
