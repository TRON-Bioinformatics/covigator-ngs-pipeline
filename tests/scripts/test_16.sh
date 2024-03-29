#!/bin/bash

##################################################################################
# Lineage only mode
##################################################################################

echo "Running CoVigator pipeline test 16"
source bin/assert.sh
output=tests/output/test15
nextflow main.nf -profile conda --name test_data \
	--output $output \
	--fasta $(pwd)/tests/test_data/test_data.fasta \
    	--lineage_mode

test -s $output/test_data.assembly.pangolin.csv || { echo "Missing pangolin output file (lineage mode with vcf input)!"; exit 1; }
assert_eq $(wc -l $output/test_data.assembly.pangolin.csv | cut -d ' ' -f1) 2 "Wrong number of pangolin results"
