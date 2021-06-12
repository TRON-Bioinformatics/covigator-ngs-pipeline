#!/bin/bash


bam=$1


samtools mpileup -aa -A -d 600000 -B -Q 0 $bam | ivar variants -p ${bam%.*}.ivar -q 20 -t 0.03 -r ../test_data/Sars_cov_2.ASM985889v3.dna.toplevel.fa -g ../test_data/Sars_cov_2.ASM985889v3.101.gff3

