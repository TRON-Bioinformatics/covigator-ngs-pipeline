#!/bin/bash

##################################################################################
# Lineage only mode
##################################################################################

echo "Running CoVigator pipeline test 15"
source bin/assert.sh
output=tests/output/test15
nextflow main.nf -profile conda --name test_data \
	--output $output \
	--vcf $(pwd)/tests/test_data/test_data.assembly.lineage_mode.vcf.gz \
        --lineage_mode \
        --skip_normalization

test -s $output/test_data.input.fasta || { echo "Missing VCF2FASTA fasta file (lineage mode with vcf input)!"; exit 1; }
test -s $output/test_data.input.pangolin.csv || { echo "Missing pangolin output file (lineage mode with vcf input)!"; exit 1; }
assert_eq $(wc -l $output/test_data.input.pangolin.csv | cut -d ' ' -f1) 2 "Wrong number of pangolin results"



    
