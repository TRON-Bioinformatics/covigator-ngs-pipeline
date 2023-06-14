process NANOCALLER {
	cpus params.cpus
	memory params.memory
	publishDir "${params.outDir}/${sampleId}/variant_calls/nanocaller", mode: "copy", pattern: "*vcf*"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::nanocaller=3.0.1" : null)

	input:
		tuple val(sampleId), path(bam), path(bai)

	output:
		tuple val(sampleId), path(bam), path(bai), path("${sampleId}.vcf.gz"), val("nanocaller"), emit: vcf
		path("${sampleId}*")
	
	"""
	NanoCaller \
	--bam ${bam} \
	--ref ${params.reference} \
	--preset ont \
	--maxcov 1000 \
	--mincov 10 \
	--cpu ${params.cpus} \
	--min_allele_freq 0.05 \
	--ins_threshold 0.05 \
	--del_threshold 0.05 \
	--sample ${sampleId} \
	--prefix ${sampleId}
	"""
}

process CLAIR3 {
	cpus params.cpus
	memory params.memory
	publishDir "${params.outDir}/${sampleId}/variant_calls", mode: "copy", pattern: "clair3"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::clair3=1.0.0" : null)

	input:
		tuple val(sampleId), path(bam), path(bai)

	output:
		tuple val(sampleId), path(bam), path(bai), path("clair3/merge_output.vcf.gz"), val("clair3"), emit: vcf
		path("clair3")
	
	"""
	run_clair3.sh \
	--bam_fn=${bam} \
	--ref_fn=${params.reference} \
	--threads=${params.cpus} \
	--platform="ont" \
	--model_path=${params.clair3model} \
	--include_all_ctgs \
	--min_coverage=10 \
	--snp_min_af=0.05 \
	--indel_min_af=0.05 \
	--output=clair3 \
	--chunk_size=5000
	
	# Quick fix for RefCalls from Clair3 which mess up Vafator run_clair3
	# Remove RefCalls from VCF
	mv clair3/merge_output.vcf.gz clair3/temp.vcf.gz
	zcat clair3/temp.vcf.gz | grep -v RefCall | gzip > clair3/merge_output.vcf.gz
	"""
}
