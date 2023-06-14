process OPTIMUS_NO_PRIME {
	cpus 2
	memory params.memory
	publishDir "${params.outDir}/${sampleId}/optimus_no_prime", mode: "copy", pattern: "*bed*"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::bedtools=2.30.0 bioconda::samtools=1.16.1" : null)

	input:
		tuple val(sampleId), path(bam), path(bai)

	output:
		tuple val(sampleId), path(bam), path(bai), env(primerset), emit: prediction
		path("${sampleId}_*")
	
	"""
	${moduleDir}/optimus_no_prime.sh ${params.primerDir} ${sampleId} ${bam}
	primerset=\$(python ${moduleDir}/optimus_no_prime.py ./ ${sampleId})
	"""
}

process PRIMER_SOFTCLIP {
	cpus 2
	memory params.memory
	publishDir "${params.outDir}/${sampleId}", mode: "copy"
	tag "${sampleId}"

	conda (params.enable_conda ? "bioconda::samtools=1.16.1" : null)

	input:
		tuple val(sampleId), path(bam), path(bai), val(bed)

	output:
		tuple val(sampleId), path("*_noprime.bam"), path("*_noprime.bam.bai"), emit: bam
		path("${sampleId}_*"), optional: true

	script:
	if( bed == 'NA' ){
		"""
		mv ${bam} ${sampleId}_noprime.bam
		mv ${bai} ${sampleId}_noprime.bam.bai
		"""
	} else {
		"""
		samtools ampliconclip \
			-b ${params.primerDir}/${bed} \
			-f ${sampleId}_ampliconclip_stats.txt \
			--both-ends \
			--strand \
			--tolerance 5 \
			${bam} \
		| samtools sort \
			-o ${sampleId}_noprime.bam \
			-T ${sampleId} \
			-

		samtools index ${sampleId}_noprime.bam
		"""
	}
}

