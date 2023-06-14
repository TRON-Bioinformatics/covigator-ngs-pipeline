process NANOFILT {
	cpus 3
	memory params.memory
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::nanofilt=2.8.0" : null)
	
	input:
		tuple val(sampleId), path(inputFq)

	output:
		tuple val(sampleId), path("${sampleId}_qctrim.fq.gz"), emit: fq
	
	"""
	gunzip -c ${inputFq} | NanoFilt \
	--quality 10 \
	--length 300 \
	--logfile ${sampleId}_qctrim.log \
	| gzip > ${sampleId}_qctrim.fq.gz
	"""
}

process PORECHOP {
	cpus params.cpus
	memory params.memory
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::porechop=0.2.4" : null)

	input:
		tuple val(sampleId), path(inputFq)

	output:
		tuple val(sampleId), path("${sampleId}_noAdapter.fq.gz"), emit: fq
		path("${sampleId}_porechop.log")
	
	"""
	porechop \
		--input ${inputFq} \
		--output ${sampleId}_noAdapter.fq.gz \
		--threads ${params.cpus} \
		--format "fastq" \
		--verbosity 1 > ${sampleId}_porechop.log
	"""
}

process CHOPPER {
	cpus params.cpus+1
	memory params.memory
	publishDir "${params.outDir}/${sampleId}", mode: "copy", pattern: "*.log"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::chopper=0.5.0" : null)

	input:
		tuple val(sampleId), path(inputFq)

	output:
		tuple val(sampleId), path("${sampleId}_qctrim.fq.gz"), emit: fq
		path("${sampleId}_chopper.log")
	
	"""
	cat ${inputFq} \
	| chopper \
		--quality 10 \
		--minlength 300 \
		--threads ${params.cpus} \
		2> ${sampleId}_chopper.log \
	| gzip > ${sampleId}_qctrim.fq.gz
	"""
}

process PORECHOP_ABI {
	cpus params.cpus
	memory params.memory
	publishDir "${params.outDir}/${sampleId}", mode: "copy", pattern: "*.log"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::porechop_abi=0.5.0" : null)

	input:
		tuple val(sampleId), path(inputFq)

	output:
		tuple val(sampleId), path("${sampleId}_noAdapter.fq.gz"), emit: fq
		path("${sampleId}_porechop_abi.log")
	
	"""
	porechop_abi \
		--ab_initio \
		--input ${inputFq} \
		--output ${sampleId}_noAdapter.fq.gz \
		--threads ${params.cpus} \
		--format "fastq.gz" \
		--verbosity 1 > "${sampleId}_porechop_abi.log"
	"""
}

process MINIMAP2 {
	cpus params.cpus
	memory params.memory
	publishDir "${params.outDir}/${sampleId}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::minimap2=2.24 bioconda::samtools=1.16.1" : null)

	input:
		tuple val(sampleId), path(fq)

	output:
		tuple val(sampleId), path("${sampleId}.bam"), path("${sampleId}.bam.bai"), emit: bam

	"""
	cores=\$((${params.cpus} - 2))

	minimap2 \
		-x map-ont \
		-a ${params.sarscov2_reference} \
		-t \${cores} \
		${fq} \
	| samtools view \
		-b \
		- \
	| samtools sort \
		-o ${sampleId}.bam \
		--threads 1 \
		-T ${sampleId} \
		-

	samtools index ${sampleId}.bam
	"""
}