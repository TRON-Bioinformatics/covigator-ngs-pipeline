params.memory = "3g"
params.cpus = 1
params.output = "."


process readTrimmingPairedEnd {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"

    input:
        tuple val(name), file(fastq1), file(fastq2)

    output:
        tuple val(name), file("${fastq1.baseName}.trimmed.fq.gz"), file("${fastq2.baseName}.trimmed.fq.gz")
        file("${name}.fastp_stats.json")
        file("${name}.fastp_stats.html")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp \
    --in1 ${fastq1} \
    --in2 ${fastq2} \
    --out1 ${fastq1.baseName}.trimmed.fq.gz \
    --out2 ${fastq2.baseName}.trimmed.fq.gz \
    --json ${name}.fastp_stats.json \
    --html ${name}.fastp_stats.html
    """
}

process readTrimmingSingleEnd {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"

    input:
        tuple val(name), file(fastq1)

    output:
        tuple val(name), file("${fastq1.baseName}.trimmed.fq.gz")
        file("${name}.fastp_stats.json")
        file("${name}.fastp_stats.html")

    """
    # --input_files needs to be forced, otherwise it is inherited from profile in tests
    fastp \
    --in1 ${fastq1} \
    --out1 ${fastq1.baseName}.trimmed.fq.gz \
    --json ${name}.fastp_stats.json \
    --html ${name}.fastp_stats.html
    """
}