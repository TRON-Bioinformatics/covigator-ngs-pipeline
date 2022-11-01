params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false


process BAM_PREPROCESSING {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0 bioconda::sambamba=0.8.2" : null)

    input:
        tuple val(name), file(bam)
        val(reference)

    output:
        tuple val(name), file("${name}.preprocessed.bam"), file("${name}.preprocessed.bai"), emit: preprocessed_bams

    """
    mkdir tmp

    gatk CleanSam \
    --java-options '-Xmx${params.memory} -Djava.io.tmpdir=./tmp' \
    --INPUT ${bam} \
    --OUTPUT /dev/stdout | \
    gatk AddOrReplaceReadGroups \
    --java-options '-Xmx${params.memory} -Djava.io.tmpdir=./tmp' \
    --VALIDATION_STRINGENCY SILENT \
    --INPUT /dev/stdin \
    --OUTPUT ${name}.prepared.bam \
    --REFERENCE_SEQUENCE ${reference} \
    --RGPU 1 \
    --RGID 1 \
    --RGSM ${name} \
    --RGLB 1 \
    --RGPL ILLUMINA

    # removes duplicates (sorted from the alignment process)
    sambamba markdup \
        -r \
        -t ${task.cpus} \
        --tmpdir=./tmp \
        ${name}.prepared.bam ${name}.preprocessed.bam

    # removes intermediate BAM files
    rm -f ${name}.prepared.bam

    # indexes the output BAM file
    sambamba index \
        -t ${task.cpus} \
        ${name}.preprocessed.bam ${name}.preprocessed.bai
    """
}

process PRIMER_TRIMMING_IVAR {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    publishDir "${params.output}", mode: "copy", pattern: "${name}.deduplication_metrics.txt"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0 bioconda::ivar=1.3.1" : null)

    input:
        tuple val(name), file(bam), file(bai)
        file(primers)

    output:
        tuple val(name), file("${bam.baseName}.trimmed.sorted.bam"), file("${bam.baseName}.trimmed.sorted.bai"), emit: trimmed_bam

    """
    ivar trim \
    -i ${bam} \
    -b ${primers} \
    -p ${bam.baseName}.trimmed

    gatk SortSam \
    --java-options '-Xmx${params.memory}  -Djava.io.tmpdir=./tmp' \
    --INPUT ${bam.baseName}.trimmed.bam \
    --OUTPUT ${bam.baseName}.trimmed.sorted.bam \
    --SORT_ORDER coordinate

    gatk BuildBamIndex --INPUT ${bam.baseName}.trimmed.sorted.bam
    """
}

process COVERAGE_ANALYSIS {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy", pattern: "${name}.coverage.tsv"
    publishDir "${params.output}", mode: "copy", pattern: "${name}.depth.tsv"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
        tuple val(name), file(bam), file(bai)

    output:
        file("${name}.coverage.tsv")
        file("${name}.depth.tsv")


    """
    samtools coverage ${bam} > ${name}.coverage.tsv
    samtools depth -s -d 0 -H ${bam} > ${name}.depth.tsv
    """
}
