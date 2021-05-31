#!/usr/bin/env nextflow

params.help= false
params.fastq1 = false
params.fastq2 = false
params.reference = false
params.output = 'output'
params.memory = "3g"
params.cpus = 1

if (params.help) {
    log.info params.help_message
    exit 0
}

if (!params.reference) {
    log.error "--reference is required"
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
    publishDir "${params.output}/${name}", mode: "copy"

    input:
        val name from params.name
    	val fastq1 from params.fastq1
    	val fastq2 from params.fastq2

    output:
	    set val("${name}"), file("${name}.bam") into bam_files

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