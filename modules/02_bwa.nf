params.memory = "3g"
params.cpus = 1


process alignmentPairedEnd {
    cpus params.cpus
    memory params.memory
    tag params.name

    input:
        tuple val(name), file(fastq1), file(fastq2)
        val(reference)

    output:
        tuple val(name), file("${name}.bam")

    """
    bwa mem -t ${task.cpus} ${reference} ${fastq1} ${fastq2} | \
    samtools view -uS - | \
    samtools sort - > ${name}.bam
    """
}

process alignmentSingleEnd {
    cpus params.cpus
    memory params.memory
    tag params.name

    input:
        tuple val(name), file(fastq1)
        val(reference)

    output:
        tuple val(name), file("${name}.bam")

    """
    bwa mem -t ${task.cpus} ${reference} ${fastq1} | \
    samtools view -uS - | \
    samtools sort - > ${name}.bam
    """
}