params.cpus = 1
params.memory = "3g"
params.output = "."
params.min_mapping_quality = 20
params.min_base_quality = 20
params.enable_conda = false


process BGZIP {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    publishDir "${params.output}", mode: "copy"

    conda (params.enable_conda ? "bioconda::samtools=1.12" : null)

    input:
    tuple val(name), file(vcf)

    output:
    tuple val(name), file("${vcf}.gz"), file("${vcf}.gz.tbi"), emit: compressed_vcfs

    script:
    """
    bgzip -c ${vcf} > ${vcf}.gz
    tabix -p vcf ${vcf}.gz
    """
}