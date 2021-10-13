params.memory = "3g"
params.cpus = 1
params.output = "."


process pangolinLineage {
    cpus params.cpus
    memory params.memory
    publishDir "${params.output}", mode: "copy"

    input:
    tuple val(name), file(fasta)

    output:
    file("lineage_report.csv") into lineage_report

    shell:
    '''
    pangolin --outdir . ${fasta}
    '''
}