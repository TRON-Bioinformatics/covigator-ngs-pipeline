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
	    set name, file("${name}.preprocessed.bam") into preprocessed_bams, preprocessed_bams2, preprocessed_bam3, preprocessed_bams4
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

process variantCallingBcfTools {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(bam) from preprocessed_bams

    output:
	    set name, file("${name}.bcftools.bcf") into bcftools_vcfs

    """
    bcftools mpileup -E -d 0 -A -f ${params.reference} -a AD ${bam} | bcftools call -mv --ploidy 1 -Ob -o ${name}.bcftools.bcf
    #bgzip ${name}.bcftools.vcf
    #tabix -p vcf ${name}.bcftools.vcf.gz
	"""
}

process variantCallingLofreq {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(bam) from preprocessed_bams2

    output:
	    set name, file("${name}.lofreq.bcf") into lofreq_vcfs

    """
    # TODO: use an inception here to avoid writing another BAM
    lofreq indelqual --dindel --ref ${params.reference} -o ${name}.lofreq.bam ${bam}
    lofreq call --min-bq 20 --min-alt-bq 20 --min-mq 20 --ref ${params.reference} --call-indels --out ${name}.lofreq.vcf ${name}.lofreq.bam
    # we need this conversion, first tabix index and then convert to BCF to set the contig in the header which is not
    # set by lofreq and bcftools complains about it
    bgzip ${name}.lofreq.vcf
    tabix -p vcf ${name}.lofreq.vcf.gz
    bcftools view -Ob -o ${name}.lofreq.bcf ${name}.lofreq.vcf.gz
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
        set name, file(vcf) from bcftools_vcfs.concat(lofreq_vcfs)

    output:
	    set name, file("${vcf.baseName}.normalized.vcf") into normalized_vcf_files

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    nextflow run tron-bioinformatics/tronflow-variant-normalization -r ${params.tronflow_variant_normalization_version} \
    --input_vcf ${vcf} \
    --input_files false \
    --output . \
    --reference ${params.reference} \
    -profile ${workflow.profile} \
    -work-dir ${workflow.workDir}

    mv ${vcf.baseName}/${vcf.baseName}.normalized.vcf ${vcf.baseName}.normalized.vcf
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
	    set name, file("${vcf.baseName}.annotated.vcf") into annotated_vcf_files

    """
     bcftools csq --fasta-ref ${params.reference} --gff-annot ${params.gff} ${vcf} -o ${vcf.baseName}.annotated.vcf
    """
}