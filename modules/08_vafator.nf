params.cpus = 1
params.memory = "3g"
params.output = "."
params.min_mapping_quality = 20
params.min_base_quality = 20
params.enable_conda = false
params.keep_intermediate = false


process VAFATOR {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::vafator=1.1.2" : null)

    input:
    tuple val(name), file(vcf), file(bam), file(bai)

    output:
    tuple val(name), file("${vcf.baseName}.vaf.vcf"), emit: annotated_vcf

    script:
    mq_param = params.min_mapping_quality ? "--mapping-quality " + params.min_mapping_quality : ""
    bq_param = params.min_base_quality ? "--base-call-quality " + params.min_base_quality : ""
    """
    vafator \
    --input-vcf ${vcf} \
    --output-vcf ${vcf.baseName}.vaf.vcf \
    --bam vafator ${bam} ${mq_param} ${bq_param}
    """
}