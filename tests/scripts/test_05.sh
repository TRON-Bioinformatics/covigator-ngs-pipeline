#!/bin/bash

##################################################################################
# FASTA input
# use --input_fastas_list
##################################################################################
echo "Running CoVigator pipeline test 5"
source bin/assert.sh
output=tests/output/test5
echo -e "test_data\t"`pwd`"/test_data/test_data.fasta\n" > tests/test_data/test_input.txt
nextflow main.nf -profile test,conda --input_fastas_list tests/test_data/test_input.txt \
	--output $output

test -s $output/test_data.assembly.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.assembly.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }

assert_eq `zcat $output/test_data.assembly.vcf.gz | grep -v '#' | wc -l` 13 "Wrong number of variants"
assert_eq `cat $output/test_data.assembly.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
