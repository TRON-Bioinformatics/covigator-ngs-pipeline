params.memory = "3g"
params.cpus = 1
params.output = "."


process PANGOLIN_LINEAGE {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::pangolin=4.1.2" : null)

    input:
        tuple val(name), val(caller), file(fasta)

    output:
        file("${name}.${caller}.pangolin.csv")

    when:
        // only runs pangolin on LoFreq and the assembly results
        caller == "lofreq" || caller == "assembly"

    shell:
    """
    mkdir tmp

    #--decompress-model
    pangolin \
    ${fasta} \
    --outfile ${name}.${caller}.pangolin.csv \
    --tempdir ./tmp \
    --threads ${params.cpus}
    """
}

process VCF2FASTA {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.14" : null)

    input:
        tuple val(name), val(caller), file(vcf)
        val(reference)

    output:
        tuple val(name), val(caller), file("${name}.${caller}.fasta")

    shell:
    """
    bcftools index ${vcf}

    # GATK results have all FILTER="."
    bcftools consensus --fasta-ref ${reference} \
    --include 'FILTER="PASS" | FILTER="."' \
    --output ${name}.${caller}.fasta \
    ${vcf}
    """
}
