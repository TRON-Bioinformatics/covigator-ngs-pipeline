process FASTQC {
	cpus 2
	memory "3 GB"
	publishDir "${params.outDir}/${sampleId}/qc", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
	
	input:
		tuple val(sampleId), path(fq)
	
	output:
		path("*fastqc*")

	"""
	fastqc \
		--threads 2 \
		${fq}
	"""
}

process NANOPLOT {
	cpus 2
	memory "3 GB"
	publishDir "${params.outDir}/${sampleId}/qc", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::nanoplot=1.40.2" : null)
	
	input:
		tuple val(sampleId), path(fq)

	output:
		path("*${sampleId}_*")
	
	"""
	NanoPlot \
		--fastq ${fq} \
		-p ${sampleId}_ \
		-t 1 \
		--tsv_stats
	"""
}
