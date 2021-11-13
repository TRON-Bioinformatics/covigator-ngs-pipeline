params.memory = "3g"
params.cpus = 1
params.output = "."


process PANGOLIN_LINEAGE {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"

    conda (params.enable_conda ? "bioconda::pangolin=3.1.16" : null)

    input:
        tuple val(name), file(fasta)

    output:
        file("lineage_report.csv")

    shell:
    """
    pangolin --outdir . ${fasta}
    """
}