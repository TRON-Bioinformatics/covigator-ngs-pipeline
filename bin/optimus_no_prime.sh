#!/bin/bash

primerDir=${1}
sampleId=${2}
bam=${3}

# Get coverage of read-ends across genome
	bedtools genomecov -5 -strand + -ibam ${bam} -dz \
		| awk '{OFS="\t"}{x=$2+1}{print $1,$2,x,$1":"$2"-"x":-",$3,"+"}' > temp.bed
	bedtools genomecov -5 -strand - -ibam ${bam} -dz \
		| awk '{OFS="\t"}{x=$2+1}{print $1,$2,x,$1":"$2"-"x":-",$3,"-"}' >> temp.bed
	sort -k1,1 -k2,2n temp.bed > ${sampleId}_5p_cov.bed
	cut -f1,2,3,5 temp.bed > ${sampleId}_5p_cov.bedGraph
	rm temp.bed

# Overlap read-end coverage with primer BED files
	primerList=($(find ${primerDir} -name "*.bed" -printf "%f\n"))
	for pset in ${primerList[@]}
	do
		bedtools map -a ${primerDir}/${pset} -b ${sampleId}_5p_cov.bed -s -c 5 -o sum -null 0 \
			| awk -F'\t' '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$(NF)}' \
			> ${sampleId}_${pset}
	done
