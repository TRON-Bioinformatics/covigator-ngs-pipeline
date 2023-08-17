params.memory = "3g"
params.cpus = 1
params.output = "."


process BWA_INDEX {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.12 bioconda::gatk4=4.2.0.0" : null)

    input:
    val(reference)

    output:
    path "bwamem2/index", emit: index
    file("${reference.baseName}.fai")
    file("index/${reference.baseName}.dict")
    
    script:
    memory = "${params.memory}".replaceAll(" ", "").toLowerCase()
    """
    mdkir bwamem2

    bwa-mem2 index ${reference} -p index/${reference.baseName}
    samtools faidx --fai-idx index/${reference.baseName}.fai ${reference}
    gatk CreateSequenceDictionary --REFERENCE ${reference} --OUTPUT index/${reference.baseName}.dict
    """
}

process SNPEFF_DATABASE {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    input:
    val(reference)
    val(gff)
    val(snpeff_organism)
    
    output:
    file("snpEff.config")

    script:
    memory = "${params.memory}".replaceAll(" ", "").toLowerCase()
    """
    mkdir ${snpeff_organism}
    echo ${snpeff_organism}.genome : ${snpeff_organism} > snpEff.config
    cp ${reference} ${snpeff_organism}/sequences.fa
    cp ${gff} ${snpeff_organism}/genes.gff
    snpEff build -gff3 -v ${snpeff_organism}
    """

}

