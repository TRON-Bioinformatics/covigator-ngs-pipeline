params.memory = "3g"
params.cpus = 1


process ALIGNMENT_PAIRED_END {
    cpus params.cpus
    memory params.memory
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

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

process ALIGNMENT_SINGLE_END {
    cpus params.cpus
    memory params.memory
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bwa=0.7.17 bioconda::samtools=1.12" : null)

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