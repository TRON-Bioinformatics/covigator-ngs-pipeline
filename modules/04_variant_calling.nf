params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false
params.min_mapping_quality = 20
params.min_base_quality = 20
params.match_score = 2
params.mismatch_score = -1
params.open_gap_score = -3
params.extend_gap_score = -0.1


process VARIANT_CALLING_BCFTOOLS {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bcftools=1.14" : null)

    input:
        tuple val(name), file(bam), file(bai)
        val(reference)

    output:
        tuple val(name), val("bcftools"), file("${name}.bcftools.bcf")

    """
    bcftools mpileup \
    --redo-BAQ \
    --max-depth 0 \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    --count-orphans \
    --fasta-ref ${reference} \
    --annotate AD ${bam} | \
    bcftools call \
    --multiallelic-caller \
    --variants-only \
     --ploidy 1 \
     --output-type b - > ${name}.bcftools.bcf
    """
}

process VARIANT_CALLING_LOFREQ {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bcftools=1.14 bioconda::lofreq=2.1.5" : null)

    input:
        tuple val(name), file(bam), file(bai)
        val(reference)

    output:
        tuple val(name), val("lofreq"), file("${name}.lofreq.bcf")

    """
    lofreq call \
    --min-bq ${params.min_base_quality} \
    --min-alt-bq ${params.min_base_quality} \
    --min-mq ${params.min_mapping_quality} \
    --ref ${reference} \
    --call-indels \
    <( lofreq indelqual --dindel --ref ${reference} ${bam} ) > ${name}.lofreq.vcf

    bgzip -c ${name}.lofreq.vcf > ${name}.lofreq.vcf.gz
    bcftools index ${name}.lofreq.vcf.gz
    bcftools view --output-type b ${name}.lofreq.vcf.gz > ${name}.lofreq.bcf
    """
}

process VARIANT_CALLING_GATK {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bcftools=1.14 bioconda::gatk4=4.2.0.0" : null)

    input:
        tuple val(name), file(bam), file(bai)
        val(reference)

    output:
        tuple val(name), val("gatk"), file("${name}.gatk.bcf")

    """
    mkdir tmp
    gatk HaplotypeCaller \
    --java-options '-Xmx${params.memory} -Djava.io.tmpdir=tmp' \
    --input $bam \
    --output ${name}.gatk.vcf \
    --reference ${reference} \
    --ploidy 1 \
    --min-base-quality-score ${params.min_base_quality} \
    --minimum-mapping-quality ${params.min_mapping_quality} \
    --annotation AlleleFraction

    bcftools view --output-type b ${name}.gatk.vcf > ${name}.gatk.bcf
    """
}


process VARIANT_CALLING_IVAR {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::samtools=1.12 bioconda::ivar=1.3.1" : null)

    input:
        tuple val(name), file(bam), file(bai)
        val(reference)
        file(gff)

    output:
        tuple val(name), file("${name}.ivar.tsv")

    """
    samtools mpileup \
    -aa \
    --count-orphans \
    --max-depth 0 \
    --redo-BAQ \
    --min-BQ ${params.min_base_quality} \
    --min-MQ ${params.min_mapping_quality} \
    ${bam} | \
    ivar variants \
    -p ${name}.ivar \
    -q ${params.min_base_quality} \
    -t 0.03 \
    -r ${reference} \
    -g ${gff}
    """
}

process IVAR2VCF {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bcftools=1.14 conda-forge::python=3.8.5 conda-forge::pandas=1.1.5 conda-forge::dataclasses=0.8 bioconda::pysam=0.17.0" : null)

    input:
        tuple val(name), file(tsv)
        val(reference)

    output:
        tuple val(name), val("ivar"), file("${name}.ivar.bcf")

    """
    ivar2vcf.py \
    --fasta ${reference} \
    --ivar ${tsv} \
    --output-vcf ${name}.ivar.vcf

    bcftools view --output-type b ${name}.ivar.vcf > ${name}.ivar.bcf
    """
}


process VARIANT_CALLING_ASSEMBLY {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::python=3.8.5 conda-forge::biopython=1.79" : null)

    input:
        tuple val(name), val(caller), file(fasta)
        val(reference)

    output:
        tuple val(name), val(caller), file("${name}.${caller}.vcf")

    """
    assembly_variant_caller.py \
    --fasta ${fasta} \
    --reference ${reference} \
    --output-vcf ${name}.${caller}.vcf \
    --match-score $params.match_score \
    --mismatch-score $params.mismatch_score \
    --open-gap-score $params.open_gap_score \
    --extend-gap-score $params.extend_gap_score
    """
}