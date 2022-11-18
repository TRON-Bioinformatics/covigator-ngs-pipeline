#!/bin/bash

##################################################################################
# primer trimmming
##################################################################################
echo "Running CoVigator pipeline test 10"
source bin/assert.sh
output=tests/output/test10
echo -e "test_data\t"`pwd`"/tests/test_data/test_data_1.fastq.gz\n" > tests/test_data/test_input.txt
nextflow main.nf -profile test,conda --input_fastqs_list tests/test_data/test_input.txt \
--library single --output $output \
--primers tests/test_data/SARS-CoV-2.primer.bed

test -s $output/test_data.lofreq.vcf.gz || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.json || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.fastp_stats.html || { echo "Missing VCF output file!"; exit 1; }
test -s $output/test_data.coverage.tsv || { echo "Missing coverage output file!"; exit 1; }
test -s $output/test_data.depth.tsv || { echo "Missing depth output file!"; exit 1; }
test -s $output/test_data.lofreq.pangolin.csv || { echo "Missing pangolin output file!"; exit 1; }

assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | wc -l` 3 "Wrong number of variants"
assert_eq `zcat $output/test_data.lofreq.vcf.gz | grep -v '#' | grep PASS | wc -l` 0 "Wrong number of variants"

assert_eq `cat $output/test_data.lofreq.pangolin.csv |  wc -l` 2 "Wrong number of pangolin results"
