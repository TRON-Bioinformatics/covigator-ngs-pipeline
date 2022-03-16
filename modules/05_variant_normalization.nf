params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false


process VARIANT_NORMALIZATION {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::gsl=2.7 bioconda::bcftools=1.12" : null)

    input:
        tuple val(name), val(caller), file(vcf)
        val(reference)

    output:
      tuple val(name), val(caller), file("${name}.${caller}.normalized.vcf")

    script:
    """
    # initial sort of the VCF
    bcftools sort ${vcf} | \

    # checks reference genome, decompose multiallelics, trim and left align indels
    bcftools norm --multiallelics -any --check-ref e --fasta-ref ${reference} \
    --old-rec-tag OLD_CLUMPED - | \

    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${name}.${caller}.normalized.vcf -
    """
}

process PHASING {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }
    tag "${name}"

    conda (params.enable_conda ? "conda-forge::python=3.8.5 conda-forge::pandas=1.1.5 bioconda::pysam=0.17.0 bioconda::gtfparse=1.2.1 bioconda::cyvcf2=0.30.14" : null)

    input:
        tuple val(name), val(caller), file(vcf)
        val(fasta)
        val(gtf)

    output:
      tuple val(name), val(caller), file("${name}.${caller}.phased.vcf")

    script:
    """
    phasing.py \
    --fasta ${fasta} \
    --gtf ${gtf} \
    --input-vcf ${vcf} \
    --output-vcf ${name}.${caller}.phased.vcf
    """
}