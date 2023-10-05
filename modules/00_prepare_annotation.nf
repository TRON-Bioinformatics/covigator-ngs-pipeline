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
        path("reference/sequences.fa"), emit: reference
        path("reference/sequences.fa.fai"), emit: fai
        path("reference/sequences.dict"), emit: gatk_dict
        path("reference/sequences.fa.0123")
        path("reference/sequences.fa.amb")
        path("reference/sequences.fa.ann")
        path("reference/sequences.fa.bwt.2bit.64")
        path("reference/sequences.fa.pac")
    
    script:
        memory = "${params.memory}".replaceAll(" ", "").toLowerCase()
        """
        mkdir -p reference
        cp ${reference} reference/sequences.fa
        bwa-mem2 index reference/sequences.fa
        samtools faidx reference/sequences.fa
        gatk CreateSequenceDictionary --REFERENCE reference/sequences.fa 
        """
}

process SNPEFF_DATABASE {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"
    tag "${name}"

    conda (params.enable_conda ? "bioconda::snpeff=5.0 bioconda::samtools=1.12" : null)

    input:
        val(reference)
        val(gff)
        val(snpeff_organism)
    
    output:
        path("snpeff/snpEff.config"), emit: snpeff_config
        path("snpeff/"), emit: snpeff_data

    script:
        memory = "${params.memory}".replaceAll(" ", "").toLowerCase()
        """
        mkdir -p snpeff/${snpeff_organism}
        echo ${snpeff_organism}.genome : ${snpeff_organism} > snpeff/snpEff.config
        cp ${reference} snpeff/${snpeff_organism}/sequences.fa
        cp ${gff} snpeff/${snpeff_organism}/genes.gff
        cd snpeff
        snpEff build -gff3 -v ${snpeff_organism} -dataDir .
        """

}

