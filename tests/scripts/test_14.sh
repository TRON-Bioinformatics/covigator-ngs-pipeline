#!/bin/bash

##################################################################################
# Genome generate
##################################################################################
echo "Running CoVigator pipeline test 14"
source bin/assert.sh
output=tests/output/test14
nextflow main.nf -profile conda --name test_data \
	--output $output \
	--reference_generate \
    --reference $(pwd)/reference/Sars_cov_2.ASM985889v3.dna.toplevel.fa  \
    --gff $(pwd)/reference/Sars_cov_2.ASM985889v3.101.gff3 \
    --snpeff_organism Sars_cov_2

# Test reference genome related output
test -s $output/reference/sequences.fa.fai || { echo "Missing fasta index file!"; exit 1; }
test -s $output/reference/sequences.dict || { echo "Missing GATK dict file!"; exit 1; }

# Test bwa index files are present
test -s $output/reference/sequences.fa.0123 || { echo "Missing bwa 0123 index file!"; exit 1; }
test -s $output/reference/sequences.fa.amb || { echo "Missing bwa amb index file!"; exit 1; }
test -s $output/reference/sequences.fa.ann || { echo "Missing bwa ann index file!"; exit 1; }
test -s $output/reference/sequences.fa.bwt.2bit.64 || { echo "Missing bwa bwt.2bit.64 index file!"; exit 1; }
test -s $output/reference/sequences.fa.pac || { echo "Missing bwa pac index file!"; exit 1; }


# Test snpEff output
test -s $output/snpeff/snpEff.config || { echo "Missing snpEff config file!"; exit 1; }
test -s $output/snpeff/Sars_cov_2/snpEffectPredictor.bin || { echo "Missing snpEff predictor bin file!"; exit 1; }
test -s $output/snpeff/Sars_cov_2/sequences.fa || { echo "Missing snpEff reference genome!"; exit 1; }
test -s $output/snpeff/Sars_cov_2/genes.gff || { echo "Missing snpEff reference annotation!"; exit 1; }
