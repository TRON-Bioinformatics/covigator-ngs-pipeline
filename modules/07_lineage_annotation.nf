params.memory = "3g"
params.cpus = 1
params.output = "."


process PANGOLIN_LINEAGE {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"

    conda (params.enable_conda ? "bioconda::pangolin=3.1.16" : null)

    input:
        tuple val(name), file(fasta)

    output:
        file("${fasta.baseName}.pangolin.csv")

    shell:
    """
    pangolin --outfile ${fasta.baseName}.pangolin.csv ${fasta}
    """
}

process VCF2FASTA {
    cpus params.cpus
    memory params.memory

    conda (params.enable_conda ? "bioconda::bcftools=1.14" : null)

    input:
        tuple val(name), file(vcf)
        val(reference)

    output:
        tuple val(name), file("${vcf.baseName}.fasta")

    shell:
    """
    bcftools index ${vcf}

    # GATK results have all FILTER="."
    bcftools consensus --fasta-ref ${reference} \
    --include 'FILTER="PASS" | FILTER="."' \
    --output ${vcf.baseName}.fasta \
    ${vcf}
    """
}
