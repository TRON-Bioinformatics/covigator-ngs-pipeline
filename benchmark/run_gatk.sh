#!/bin/bash


bam=$1

gatk HaplotypeCaller --input $bam --output ${bam%.*}.gatk.vcf --reference ../test_data/Sars_cov_2.ASM985889v3.dna.toplevel.fa --ploidy 1 --min-base-quality-score 20