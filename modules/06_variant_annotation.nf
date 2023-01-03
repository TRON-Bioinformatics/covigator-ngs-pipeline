params.memory = "3g"
params.cpus = 1
params.output = "."
params.conservation_sarscov2 = false
params.conservation_sarscov2_header = false
params.conservation_sarbecovirus = false
params.conservation_sarbecovirus_header = false
params.conservation_vertebrate = false
params.conservation_vertebrate_header = false
params.keep_intermediate = false
params.low_frequency_variant_threshold = 0.02
params.subclonal_variant_threshold = 0.5
params.lq_clonal_variant_threshold = 0.8
params.vafator_min_mapping_quality = 0
params.vafator_min_base_quality = 0


process VARIANT_ANNOTATION {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy", pattern: "*.vcf.gz*"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::snpeff=5.0 bioconda::samtools=1.12" : null)

    input:
    tuple val(name), val(caller), file(vcf)
    val(snpeff_data)
    val(snpeff_config)
    val(snpeff_organism)

    output:
    tuple val(name), val(caller),
        file("${name}.${caller}.vcf.gz"),
        file("${name}.${caller}.vcf.gz.tbi"), emit: annotated_vcfs

    script:
    memory = "${params.memory}".replaceAll(" ", "").toLowerCase()
    """
    # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
    cp ${snpeff_config} .

    snpEff eff -Xmx${memory} -dataDir ${snpeff_data} \
    -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
    ${snpeff_organism}  ${vcf} | bgzip -c > ${name}.${caller}.vcf.gz

    tabix -p vcf ${name}.${caller}.vcf.gz
    """
}

process VARIANT_VAF_ANNOTATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.14" : null)

    input:
        tuple val(name), val(caller), file(vcf)

    output:
        tuple val(name), val(caller), file("${name}.${caller}.vcf"), emit: vaf_annotated

    """
    bgzip -c ${vcf} > ${name}.vcf.gz

    tabix -p vcf ${name}.vcf.gz

    # annotates low frequency and subclonal variants
    bcftools view -Ob ${name}.vcf.gz | \
    bcftools filter \
    --exclude 'INFO/vafator_af < ${params.low_frequency_variant_threshold}' \
    --soft-filter LOW_FREQUENCY - | \
    bcftools filter \
    --exclude 'INFO/vafator_af >= ${params.low_frequency_variant_threshold} && INFO/vafator_af < ${params.subclonal_variant_threshold}' \
    --soft-filter SUBCLONAL \
    --output-type v - | \
    bcftools filter \
    --exclude 'INFO/vafator_af >= ${params.subclonal_variant_threshold} && INFO/vafator_af < ${params.lq_clonal_variant_threshold}' \
    --soft-filter LOW_QUALITY_CLONAL \
    --output-type v - > ${name}.${caller}.vcf
    """
}


process VARIANT_SARSCOV2_ANNOTATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.14" : null)

    input:
    tuple val(name), val(caller), file(vcf)

    output:
    tuple val(name), val(caller), file("${name}.${caller}.annotated_sarscov2.vcf"), emit: annotated_vcfs

    """
    bcftools annotate \
    --annotations ${params.conservation_sarscov2} \
    --header-lines ${params.conservation_sarscov2_header} \
    -c CHROM,FROM,TO,CONS_HMM_SARS_COV_2 \
    --output-type z ${vcf} | \
    bcftools annotate \
    --annotations ${params.conservation_sarbecovirus} \
    --header-lines ${params.conservation_sarbecovirus_header} \
    -c CHROM,FROM,TO,CONS_HMM_SARBECOVIRUS \
    --output-type z - | \
    bcftools annotate \
    --annotations ${params.conservation_vertebrate} \
    --header-lines ${params.conservation_vertebrate_header} \
    -c CHROM,FROM,TO,CONS_HMM_VERTEBRATE_COV \
    --output-type z - | \
    bcftools annotate \
    --annotations ${params.pfam_names} \
    --header-lines ${params.pfam_names_header} \
    -c CHROM,FROM,TO,PFAM_NAME \
    --output-type z - | \
    bcftools annotate \
    --annotations ${params.pfam_descriptions} \
    --header-lines ${params.pfam_descriptions_header} \
    -c CHROM,FROM,TO,PFAM_DESCRIPTION - > ${name}.${caller}.annotated_sarscov2.vcf

    # TODO: include this step for FASTA data
    #bcftools annotate \
    #--annotations ${params.problematic_sites} \
    #--columns FILTER \
    #--output-type b - > ${vcf.baseName}.annotated.vcf.gz
    """
}

process VAFATOR {
    cpus params.cpus
    memory params.memory
    tag "${name}"
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::vafator=1.2.5" : null)

    input:
    tuple val(name), val(caller), file(vcf), file(bam), file(bai)

    output:
    tuple val(name), val(caller), file("${name}.${caller}.vaf.vcf"), emit: annotated_vcf

    script:
    mq_param = params.vafator_min_mapping_quality != false ? "--mapping-quality " + params.vafator_min_mapping_quality : ""
    bq_param = params.vafator_min_base_quality != false ? "--base-call-quality " + params.vafator_min_base_quality : ""
    """
    vafator \
    --input-vcf ${vcf} \
    --output-vcf ${name}.${caller}.vaf.vcf \
    --bam vafator ${bam} ${mq_param} ${bq_param}
    """
}
