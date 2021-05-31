#!/usr/bin/env nextflow

params.help= false
params.fastq1 = false
params.fastq2 = false
params.reference = false
params.gff = false
params.output = false
params.memory = "3g"
params.cpus = 1

if (params.help) {
    log.info params.help_message
    exit 0
}

if (!params.output) {
    log.error "--output is required"
    exit 1
}

if (!params.reference) {
    log.error "--reference is required"
    exit 1
}

if (!params.gff) {
    log.error "--gff is required"
    exit 1
}

if (!params.name) {
    log.error "--name is required"
    exit 1
}

if (!params.fastq1) {
    log.error "--fastq1 is required"
    exit 1
}

library = "paired"
if (!params.fastq2) {
    params.fastq2 = ""
    library = "single"
}

process alignment {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
        val name from params.name
    	val fastq1 from params.fastq1
    	val fastq2 from params.fastq2

    output:
	    set val("${name}"), val("${params.output}/${name}.bam") into bam_files
	    file "${name}.bam"

    """
    echo "${name}\t${fastq1}\t${fastq2}" > bwa_input_files.txt

    nextflow run tron-bioinformatics/tronflow-bwa -r ${params.tronflow_bwa_version} \
    --input_files bwa_input_files.txt \
    --algorithm mem \
    --library ${library} \
    --output . \
    --reference ${params.reference} \
    --cpus ${task.cpus} --memory ${task.memory} \
    -profile ${workflow.profile} \
    -work-dir ${workflow.workDir}
	"""
}

process bamPreprocessing {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
        set name, bam from bam_files

    output:
	    set val("${name}"), val("${params.output}/${name}.preprocessed.bam") into preprocessed_bam_files
	    file "${name}.preprocessed.bam"
	    file "${name}.preprocessed.bai"

    """
    echo "${name}\tnormal\t${bam}" > bam_preprocessing_input_files.txt

    nextflow run tron-bioinformatics/tronflow-bam-preprocessing -r ${params.tronflow_bam_preprocessing_version} \
    --input_files bam_preprocessing_input_files.txt \
    --output . \
    --reference ${params.reference} \
    --skip_bqsr --skip_realignment --skip_metrics \
    --prepare_bam_cpus ${params.cpus} --prepare_bam_memory ${params.memory} \
    --mark_duplicates_cpus ${params.cpus} --mark_duplicates_memory ${params.memory} \
    -profile ${workflow.profile} \
    -work-dir ${workflow.workDir}

    mv ${name}/${name}.preprocessed.bam ${name}.preprocessed.bam
    mv ${name}/${name}.preprocessed.bai ${name}.preprocessed.bai
	"""
}

process variantCalling {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
        set name, bam from preprocessed_bam_files

    output:
	    set val("${name}"), val("${params.output}/${name}.vcf") into vcf_files
	    file "${name}.vcf"

    """
    bcftools mpileup -E -d 0 -A -f ${params.reference} -a AD ${bam} | bcftools call -mv --ploidy 1 -Ov -o ${name}.vcf
	"""
}

process variantNormalization {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
        set name, vcf from vcf_files

    output:
	    set val("${name}"), val("${params.output}/${name}.normalized.vcf") into normalized_vcf_files
	    file "${name}.normalized.vcf"

    """
    echo "${name}\t${vcf}" > vcf_normalization_input_files.txt

    nextflow run tron-bioinformatics/tronflow-variant-normalization -r ${params.tronflow_variant_normalization_version} \
    --input_files vcf_normalization_input_files.txt \
    --output . \
    --reference ${params.reference} \
    -profile ${workflow.profile} \
    -work-dir ${workflow.workDir}

    mv ${name}/${name}.normalized.vcf ${name}.normalized.vcf
	"""
}

process variantAnnotation {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    input:
        set name, vcf from normalized_vcf_files

    output:
	    set val("${name}"), val("${params.output}/${name}.vcf") into annotated_vcf_files
	    file "${name}.annotated.vcf"

    """
     bcftools csq --fasta-ref ${params.reference} --gff-annot ${params.gff} ${vcf} -o ${name}.annotated.vcf
    """
}