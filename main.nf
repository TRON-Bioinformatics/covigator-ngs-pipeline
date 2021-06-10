#!/usr/bin/env nextflow

params.help= false
params.fastq1 = false
params.fastq2 = false
params.name = false
params.reference = false
params.gff = false
params.output = false
params.min_mapping_quality = 20
params.min_base_quality = 20
params.low_frequency_variant_threshold = 0.2
params.subclonal_variant_threshold = 0.8
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
        nextflow run tron-bioinformatics/tronflow-bwa \
        -r ${params.tronflow_bwa_version} \
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
	    set name, file("${name}.preprocessed.bam"), file("${name}.preprocessed.bai") into preprocessed_bams,
	        preprocessed_bams2, preprocessed_bams3, preprocessed_bams4, preprocessed_bams5


    """
    # --input_files, --known_indels1 and --known_indels2 needs to be forced, otherwise it is inherited from test profile
    nextflow run tron-bioinformatics/tronflow-bam-preprocessing \
    -r ${params.tronflow_bam_preprocessing_version} \
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
        set name, file(bam), file(bai) from preprocessed_bams

    output:
	    set name, file("${name}.bcftools.bcf") into bcftools_vcfs

    """
    bcftools mpileup \
    --redo-BAQ \
    --max-depth 0 \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    --count-orphans \
    --fasta-ref ${params.reference} \
    --annotate AD ${bam} | \
    bcftools call \
    --multiallelic-caller \
    --variants-only \
     --ploidy 1 \
     --output-type b \
     --output ${name}.bcftools.bcf
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
        set name, file(bam), file(bai) from preprocessed_bams2

    output:
	    set name, file("${name}.lofreq.vcf") into lofreq_vcfs

    """
    lofreq call \
    --min-bq ${params.min_base_quality} \
    --min-alt-bq ${params.min_base_quality} \
    --min-mq ${params.min_mapping_quality} \
    --ref ${params.reference} \
    --call-indels \
    <( lofreq indelqual --dindel --ref ${params.reference} ${bam} ) | \
    bgzip -c > ${name}.lofreq.vcf.gz

    tabix -p vcf ${name}.lofreq.vcf.gz

    bcftools view -Ob -o ${name}.lofreq.bcf ${name}.lofreq.vcf.gz

    # annotates low frequency and subclonal variants
    bcftools filter \
    --exclude 'INFO/AF < ${params.low_frequency_variant_threshold}' \
    --soft-filter LOW_FREQUENCY \
    ${name}.lofreq.bcf > ${name}.lofreq2.vcf

    bcftools filter \
    --exclude 'INFO/AF >= ${params.low_frequency_variant_threshold} && INFO/AF < ${params.subclonal_variant_threshold}' \
    --soft-filter SUBCLONAL \
    ${name}.lofreq2.vcf > ${name}.lofreq.vcf
	"""
}

process variantCallingGatk {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(bam), file(bai) from preprocessed_bams3

    output:
	    set name, file("${name}.gatk.vcf") into gatk_vcfs

    """
    gatk HaplotypeCaller \
    --input $bam \
    --output ${name}.gatk.vcf \
    --reference ${params.reference} \
    --ploidy 1 \
    --min-base-quality-score ${params.min_base_quality} \
    --minimum-mapping-quality ${params.min_mapping_quality} \
    --annotation AlleleFraction
	"""
}

process variantCallingIvar {
    cpus params.cpus
    memory params.memory
    tag params.name
    publishDir "${params.output}/${params.name}", mode: "copy"

    input:
        set name, file(bam), file(bai) from preprocessed_bams4

    output:
	    file("${name}.ivar.tsv")

    """
    samtools mpileup \
    -aa \
    --count-orphans \
    --max-depth 0 \
    --redo-BAQ \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    ${bam} | \
    ivar variants \
    -p ${name}.ivar \
    -q ${params.min_base_quality} \
    -t 0.03 \
    -r ${params.reference} \
    -g ${params.gff}
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
        set name, file(vcf) from bcftools_vcfs.concat(lofreq_vcfs).concat(gatk_vcfs)

    output:
	    set name, file("${vcf.baseName}.normalized.vcf") into normalized_vcf_files

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    nextflow run tron-bioinformatics/tronflow-variant-normalization \
    -r ${params.tronflow_variant_normalization_version} \
    --input_vcf ${vcf} \
    --input_files false \
    --output . \
    --reference ${params.reference} \
    -profile ${workflow.profile} \
    -work-dir ${workflow.workDir}

    mv ${vcf.baseName}/${vcf.baseName}.normalized.vcf ${vcf.baseName}.normalized.vcf
	"""
}

process phasing {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}/${params.name}", mode: "copy"
    }

    input:
        set name, file(vcf), file(bam), file(bai) from normalized_vcf_files.combine(preprocessed_bams5, by:0)

    output:
	    set name, file("${vcf.baseName}.phased.vcf") into phased_variants

    """
    whatshap polyphase \
    --ploidy 1 \
    --indels \
    --mapping-quality ${params.min_mapping_quality} \
    --output ${vcf.baseName}.phased.vcf \
    ${vcf} \
    ${bam}
    """
}

process variantAnnotation {
    cpus params.cpus
    memory params.memory
    tag params.name
    publishDir "${params.output}/${params.name}", mode: "copy"

    input:
        set name, file(vcf) from phased_variants

    output:
	    file("${vcf.baseName}.annotated.vcf.gz")
	    file("${vcf.baseName}.annotated.vcf.gz.tbi")

    """
     bcftools csq --fasta-ref ${params.reference} --gff-annot ${params.gff} ${vcf} | \
     bgzip -c > ${vcf.baseName}.annotated.vcf.gz

     tabix -p vcf ${vcf.baseName}.annotated.vcf.gz
    """
}
