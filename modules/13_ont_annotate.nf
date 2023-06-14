process VAFATOR {
    cpus 1
    memory "3 GB"
	tag "${sampleId}"
	publishDir "${params.outDir}/${sampleId}/variant_calls/${caller}", mode: "copy"	

	conda (params.enable_conda ? "bioconda::vafator=1.2.5" : null)

	input:
	tuple val(sampleId), path(bam), path(bai), path(vcf), val(caller)
	
	output:
	tuple val(sampleId), path("${sampleId}_${caller}_vaf.vcf"), val(caller), emit: annotated_vcf

	script:
	mq_param = params.vafator_min_mapping_quality != false ? "--mapping-quality " + params.vafator_min_mapping_quality : ""
	bq_param = params.vafator_min_base_quality != false ? "--base-call-quality " + params.vafator_min_base_quality : ""
	"""
	vafator \
	--input-vcf ${vcf} \
	--output-vcf ${sampleId}_${caller}_vaf.vcf \
	--bam vafator ${bam} ${mq_param} ${bq_param}
	"""
}

process SNPEFF {
    cpus 1
    memory params.memory
	publishDir "${params.outDir}/${sampleId}/variant_calls/${caller}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::snpeff=5.0 bioconda::samtools=1.16.1" : null)

	input:
	tuple val(sampleId), path(vcf), val(caller)
	val(snpeff_data)
	val(snpeff_config)
	val(snpeff_organism)

	output:
	tuple val(sampleId), path("${sampleId}_${caller}_snpeff.vcf.gz"), emit: annotated_vcf
	path("${sampleId}_${caller}_snpeff.vcf.gz.tbi")

	"""
	# for some reason the snpEff.config file needs to be in the folder where snpeff runs...
	cp ${snpeff_config} .

	snpEff eff -Xmx${params.memory} -dataDir ${snpeff_data} \
	-noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
	${snpeff_organism} ${vcf} | bgzip -c > ${sampleId}_${caller}_snpeff.vcf.gz

	tabix -p vcf ${sampleId}_${caller}_snpeff.vcf.gz
	"""
}


process VARIANT_VAF_ANNOTATION {
	cpus 1
	memory "3 GB"
	publishDir "${params.outDir}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.14" : null)

	input:
		tuple val(sampleId), path(vcf)

	output:
		tuple val(sampleId), path("${sampleId}.lowfreq_subclonal.vcf"), emit: vaf_annotated

	"""
	tabix -p vcf ${vcf}

	# annotates low frequency and subclonal variants
	bcftools view -Ob ${vcf} | \
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
	--output-type v - > ${sampleId}.lowfreq_subclonal.vcf
	"""
}

/**
Add SARS-CoV-2 specific annotations: conservation, Pfam domains and problematic sites.
Also, according to problematic sites described in DeMaio et al. (2020) we filter out any variants at the beginning and
end of the genome.
*/
process VARIANT_SARSCOV2_ANNOTATION {
	cpus 1
	memory "3 GB"
	publishDir "${params.output}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.14" : null)

	input:
	tuple val(sampleId), path(vcf)

	output:
	tuple val(sampleId), path("${sampleId}.annotated_sarscov2.vcf"), emit: annotated_vcfs

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
	-c CHROM,FROM,TO,PFAM_DESCRIPTION - | \
	bcftools filter \
	--exclude 'POS <= 55 | POS >= 29804' \
	--output-type z - > annotated_sarscov2.vcf.gz

	tabix -p vcf annotated_sarscov2.vcf.gz

	bcftools annotate \
	--annotations ${params.problematic_sites} \
	--columns INFO/problematic:=FILTER annotated_sarscov2.vcf.gz > ${sampleId}.annotated_sarscov2.vcf
	"""
}

process PANGOLIN {
    cpus 1
    memory "3 GB"
	publishDir "${params.outDir}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::pangolin=4.1.2" : null)

	input:
	tuple val(sampleId), path(fasta)

	output:
	tuple val(name), path("*pangolin.csv")

    """
    mkdir tmp
    #--decompress-model
    pangolin \
    ${fasta} \
    --outfile ${sampleId}.pangolin.csv \
    --tempdir ./tmp \
    --threads 1
    """
}