#!/bin/bash


bam=$1

lofreq indelqual --dindel --ref ../test_data/Sars_cov_2.ASM985889v3.dna.toplevel.fa -o ${bam%.*}.lofreq.bam $bam
lofreq call --min-bq 20 --min-alt-bq 20 --min-mq 20 --ref ../test_data/Sars_cov_2.ASM985889v3.dna.toplevel.fa --call-indels --out ${bam%.*}.lofreq.vcf ${bam%.*}.lofreq.bam