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


process variantAnnotation {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"

    input:
        tuple val(name), file(vcf)

    output:
        file("${vcf.baseName}.annotated.vcf.gz")
        file("${vcf.baseName}.annotated.vcf.gz.tbi")

    """
    # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
    cp ${params.snpeff_config} .

    snpEff eff -dataDir ${params.snpeff_data} \
    -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
    ${params.snpeff_organism}  ${vcf} | \
    bgzip -c > ${vcf.baseName}.annotated.vcf.gz

    tabix -p vcf ${vcf.baseName}.annotated.vcf.gz
    """
}


process variantSarsCov2Annotation {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"

    input:
        tuple val(name), file(vcf)

    output:
        file("${vcf.baseName}.annotated.vcf.gz")
        file("${vcf.baseName}.annotated.vcf.gz.tbi")

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
    -c CHROM,FROM,TO,PFAM_DESCRIPTION \
    --output-type z - > ${vcf.baseName}.annotated.vcf.gz

    # TODO: include this step for GISAID data
    #bcftools annotate \
    #--annotations ${params.problematic_sites} \
    #--columns FILTER \
    #--output-type b - > ${vcf.baseName}.annotated.vcf.gz

    tabix -p vcf ${vcf.baseName}.annotated.vcf.gz
    """
}