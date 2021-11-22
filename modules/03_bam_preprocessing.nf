params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false


process BAM_PREPROCESSING {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    publishDir "${params.output}", mode: "copy", pattern: "${name}.deduplication_metrics.txt"

    conda (params.enable_conda ? "bioconda::gatk4=4.2.0.0" : null)

    input:
        tuple val(name), file(bam)
        val(reference)

    output:
        tuple val(name), file("${name}.preprocessed.bam"), file("${name}.preprocessed.bai"), emit: preprocessed_bam
        file("${name}.deduplication_metrics.txt")

    """
    gatk CleanSam \
    --java-options '-Xmx${params.memory} -Djava.io.tmpdir=tmp' \
    --INPUT ${bam} \
    --OUTPUT /dev/stdout | \
    gatk AddOrReplaceReadGroups \
    --java-options '-Xmx${params.memory} -Djava.io.tmpdir=tmp' \
    --VALIDATION_STRINGENCY SILENT \
    --INPUT /dev/stdin \
    --OUTPUT ${bam.baseName}.prepared.bam \
    --REFERENCE_SEQUENCE ${reference} \
    --RGPU 1 \
    --RGID 1 \
    --RGSM ${name} \
    --RGLB 1 \
    --RGPL ILLUMINA \
    --SORT_ORDER queryname

    gatk MarkDuplicates \
    --java-options '-Xmx${params.memory}  -Djava.io.tmpdir=tmp' \
    --INPUT ${bam.baseName}.prepared.bam \
    --METRICS_FILE ${name}.deduplication_metrics.txt \
    --OUTPUT ${bam.baseName}.dedup.bam \
    --REMOVE_DUPLICATES true

    gatk SortSam \
    --java-options '-Xmx${params.memory}  -Djava.io.tmpdir=tmp' \
    --INPUT ${bam.baseName}.dedup.bam \
    --OUTPUT ${bam.baseName}.preprocessed.bam \
    --SORT_ORDER coordinate

    gatk BuildBamIndex --INPUT ${bam.baseName}.preprocessed.bam
    """
}

process COVERAGE_ANALYSIS {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy", pattern: "${name}.coverage.tsv"
    publishDir "${params.output}", mode: "copy", pattern: "${name}.depth.tsv"

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