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

    conda (params.enable_conda ? "bioconda::vt=0.57721 bioconda::bcftools=1.12" : null)

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

    # decompose complex variants
    vt decompose_blocksub -a -p - | \

    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${name}.${caller}.normalized.vcf -
    """
}