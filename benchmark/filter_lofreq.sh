#!/bin/bash

input_vcf=$1


bcftools filter --exclude "INFO/AF < 0.2" $input_vcf > ${input_vcf%.*}.filter02.vcf
bcftools filter --exclude "INFO/AF < 0.1" $input_vcf > ${input_vcf%.*}.filter01.vcf
