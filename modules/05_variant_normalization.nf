params.memory = "3g"
params.cpus = 1
params.output = "."
params.keep_intermediate = false


process variantNormalization {
    cpus params.cpus
    memory params.memory
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    input:
        tuple val(name), file(vcf)
        val(reference)

    output:
      tuple val(name), file("${vcf.baseName}.normalized.vcf")

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
    bcftools norm --rm-dup exact -o ${vcf.baseName}.normalized.vcf -
    """
}