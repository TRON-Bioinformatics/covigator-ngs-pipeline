params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false


process bamPreprocessing {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    publishDir "${params.output}", mode: "copy", pattern: "${name}.deduplication_metrics.txt"
    publishDir "${params.output}", mode: "copy", pattern: "${name}.coverage.tsv"
    publishDir "${params.output}", mode: "copy", pattern: "${name}.depth.tsv"

    input:
        tuple val(name), file(bam)
        val(reference)

    output:
        tuple val(name), file("${name}.preprocessed.bam"), file("${name}.preprocessed.bai")
        file("${name}.deduplication_metrics.txt")
        file("${name}.coverage.tsv")
        file("${name}.depth.tsv")


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
    --OUTPUT ${bam.baseName}.dedup.sorted.bam \
    --SORT_ORDER coordinate

    gatk BuildBamIndex --INPUT ${bam.baseName}.dedup.sorted.bam

    gatk3 -Xmx${params.memory} -Djava.io.tmpdir=tmp -T RealignerTargetCreator \
    --input_file ${bam.baseName}.dedup.sorted.bam \
    --out ${bam.baseName}.RA.intervals \
    --reference_sequence ${reference}

    gatk3 -Xmx${params.memory} -Djava.io.tmpdir=tmp -T IndelRealigner \
    --input_file ${bam.baseName}.dedup.sorted.bam \
    --out ${name}.preprocessed.bam \
    --reference_sequence ${reference} \
    --targetIntervals ${bam.baseName}.RA.intervals \
    --consensusDeterminationModel USE_SW \
    --LODThresholdForCleaning 0.4 \
    --maxReadsInMemory 600000

    samtools coverage ${name}.preprocessed.bam > ${name}.coverage.tsv

    samtools depth -s -d 0 -H ${name}.preprocessed.bam > ${name}.depth.tsv
    """
}