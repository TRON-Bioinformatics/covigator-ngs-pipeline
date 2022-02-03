params.memory = "3g"
params.cpus = 1
params.output = "."
params.snpeff_data = false
params.snpeff_config = false
params.snpeff_organism = false
params.conservation_sarscov2 = false
params.conservation_sarscov2_header = false
params.conservation_sarbecovirus = false
params.conservation_sarbecovirus_header = false
params.conservation_vertebrate = false
params.conservation_vertebrate_header = false
params.keep_intermediate = false
params.low_frequency_variant_threshold = 0.2
params.subclonal_variant_threshold = 0.8


process VARIANT_ANNOTATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::snpeff=5.0" : null)

    input:
    tuple val(name), val(caller), file(vcf)

    output:
    tuple val(name), val(caller), file("${name}.${caller}.annotated.vcf"), emit: annotated_vcfs

    """
    # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
    cp ${params.snpeff_config} .

    snpEff eff -dataDir ${params.snpeff_data} \
    -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
    ${params.snpeff_organism}  ${vcf} > ${name}.${caller}.annotated.vcf
    """
}

process VARIANT_VAF_ANNOTATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::bcftools=1.14" : null)

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
    --output-type v - > ${name}.${caller}.vcf
    """
}


process VARIANT_SARSCOV2_ANNOTATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    conda (params.enable_conda ? "bioconda::snpeff=5.0 bioconda::bcftools=1.12" : null)

    input:
    tuple val(name), val(caller), file(vcf)

    output:
    tuple val(name), val(caller), file("${name}.${caller}.annotated.vcf"), emit: annotated_vcfs

    """
    # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
    cp ${params.snpeff_config} .

    snpEff eff -dataDir ${params.snpeff_data} \
    -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
    Sars_cov_2.ASM985889v3.101  ${vcf} | \
    bgzip -c | \
    bcftools annotate \
    --annotations ${params.conservation_sarscov2} \
    --header-lines ${params.conservation_sarscov2_header} \
    -c CHROM,FROM,TO,CONS_HMM_SARS_COV_2 \
    --output-type z - | \
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
    -c CHROM,FROM,TO,PFAM_DESCRIPTION - > ${name}.${caller}.annotated.vcf

    # TODO: include this step for GISAID data
    #bcftools annotate \
    #--annotations ${params.problematic_sites} \
    #--columns FILTER \
    #--output-type b - > ${vcf.baseName}.annotated.vcf.gz
    """
}
