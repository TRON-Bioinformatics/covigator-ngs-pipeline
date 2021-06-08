#!/usr/bin/env nextflow

params.help= false
params.fastq1 = false
params.fastq2 = false
params.name = false
params.reference = false
params.gff = false
params.output = false
params.memory = "3g"
params.cpus = 1
params.keep_intermediate = false

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
    library = "single"
}

if (library == "paired") {
    process alignmentPairedEnd {
        cpus params.cpus
        memory params.memory
        tag params.name

        input:
            val name from params.name
            file fastq1 from file(params.fastq1)
            file fastq2 from file(params.fastq2)

        output:
            set name, file("${name}.bam") into bam_files

        """
        # --input_files needs to be forced, otherwise it is inherited from profile in tests
        nextflow run tron-bioinformatics/tronflow-bwa -r ${params.tronflow_bwa_version} \
        --input_name ${name} \
        --input_fastq1 ${fastq1} \
        --input_fastq2 ${fastq2} \
        --input_files false \
        --algorithm mem \
        --library ${library} \
        --output . \
        --reference ${params.reference} \
        --cpus ${task.cpus} --memory ${task.memory} \
        -profile ${workflow.profile} \
        -work-dir ${workflow.workDir}
        """
    }
}
else {
    process alignmentSingleEnd {
        cpus params.cpus
        memory params.memory
        tag params.name

        input:
            val name from params.name
            file fastq1 from file(params.fastq1)

        output:
            set name, file("${name}.bam") into bam_files

        """
        # --input_files needs to be forced, otherwise it is inherited from profile in tests
        nextflow run tron-bioinformatics/tronflow-bwa -r ${params.tronflow_bwa_version} \
        --input_name ${name} \
        --input_fastq1 ${fastq1} \
        --input_files false \
        --algorithm mem \
        --library ${library} \
        --output . \
        --reference ${params.reference} \
        --cpus ${task.cpus} --memory ${task.memory} \
        -profile ${workflow.profile} \
        -work-dir ${workflow.workDir}
        """
    }
}

process bamPreprocessing {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(bam) from bam_files

    output:
	    set name, file("${name}.preprocessed.bam") into preprocessed_bam_files
	    file "${name}.preprocessed.bai"

    """
    # --input_files, --known_indels1 and --known_indels2 needs to be forced, otherwise it is inherited from test profile
    nextflow run tron-bioinformatics/tronflow-bam-preprocessing -r ${params.tronflow_bam_preprocessing_version} \
    --input_bam ${bam} \
    --input_files false \
    --output . \
    --reference ${params.reference} \
    --skip_bqsr --skip_metrics \
    --known_indels1 false --known_indels2 false \
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
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(bam) from preprocessed_bam_files

    output:
	    set name, file("${name}.vcf") into vcf_files

    """
    bcftools mpileup -E -d 0 -A -f ${params.reference} -a AD ${bam} | bcftools call -mv --ploidy 1 -Ov -o ${name}.vcf
	"""
}

process variantNormalization {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(vcf) from vcf_files

    output:
	    set name, file("${name}.normalized.vcf") into normalized_vcf_files

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    nextflow run tron-bioinformatics/tronflow-variant-normalization -r ${params.tronflow_variant_normalization_version} \
    --input_vcf ${vcf} \
    --input_files false \
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
    tag params.name
    publishDir "${params.output}/${params.name}", mode: "copy"

    input:
        set name, file(vcf) from normalized_vcf_files

    output:
	    set name, file("${name}.annotated.vcf") into annotated_vcf_files

    """
     bcftools csq --fasta-ref ${params.reference} --gff-annot ${params.gff} ${vcf} -o ${name}.annotated.vcf
    """
}