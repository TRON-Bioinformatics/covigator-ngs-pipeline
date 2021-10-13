#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { readTrimmingPairedEnd; readTrimmingSingleEnd } from './modules/01_fastp'
include { alignmentPairedEnd; alignmentSingleEnd } from './modules/02_bwa'
include { bamPreprocessing } from './modules/03_bam_preprocessing'
include { variantCallingBcfTools; variantCallingLofreq ; variantCallingGatk ; variantCallingIvar ; assemblyVariantCaller } from './modules/04_variant_calling'
include { variantNormalization } from './modules/05_variant_normalization'
include { variantAnnotation; variantSarsCov2Annotation } from './modules/06_variant_annotation'


params.help= false
params.initialize = false
if (params.initialize) {
    params.fastq1 = "$baseDir/test_data/ERR4145453_1.fastq.gz"
    params.skip_bcftools = true
    params.skip_ivar = true
    params.skip_gatk = true
    params.name = "init"
}
else {
    params.fastq1 = false
    params.skip_ivar = false
    params.skip_bcftools = false
    params.skip_gatk = false
    params.name = false
}

params.skip_lofreq = false
params.fasta = false
params.fastq2 = false
params.reference = false
params.gff = false
params.output = "."
params.min_mapping_quality = 20
params.min_base_quality = 20
params.low_frequency_variant_threshold = 0.2
params.subclonal_variant_threshold = 0.8
params.memory = "3g"
params.cpus = 1
params.keep_intermediate = false
params.match_score = 2
params.mismatch_score = -1
params.open_gap_score = -3
params.extend_gap_score = -0.1
params.chromosome = "MN908947.3"
params.skip_sarscov2_annotations = false
params.library = false
params.input_fastqs_list = false
params.input_fastas_list = false

if (params.help) {
    log.info params.help_message
    exit 0
}
if (!params.output) {
    log.error "--output is required"
    exit 1
}
if (!params.reference) {
    log.error "--reference is required"
    exit 1
}
if (params.fastq1 && params.fasta) {
    log.error "provide only --fastq1 or --fasta"
    exit 1
}
if (params.input_fastqs_list && params.input_fastas_list) {
    log.error "provide only --input_fastqs_list or --input_fastas_list"
    exit 1
}

input_fastqs = false
input_fastas = false
if (params.input_fastqs_list || params.fastq1) {

    if (!params.gff) {
        log.error "--gff is required"
        exit 1
    }
    else {
        gff = file(params.gff)
    }

    // if independent FASTQ files are provided the value of library is overridden
    if (params.fastq1 && params.fastq2) {
        params.library = "single"
    }
    else if (params.fastq1 && params.fastq2) {
        params.library = "paired"
    }
    else if (params.input_fastqs_list && !params.library) {
        log.error "--library paired|single is required when --input_fastqs_list is provided"
        exit 1
    }

    if (params.input_fastqs_list) {
        if (params.library == "paired") {
            Channel
                .fromPath(params.input_fastqs_list)
                .splitCsv(header: ['name', 'fastq1', 'fastq2'], sep: "\t")
                .map{ row-> tuple(row.name, file(row.fastq1), file(row.fastq2)) }
                .set { input_fastqs }
        }
        else {
            Channel
                .fromPath(params.input_fastqs_list)
                .splitCsv(header: ['name', 'fastq'], sep: "\t")
                .map{ row-> tuple(row.name, file(row.fastq)) }
                .set { input_fastqs }
        }
    }
    else {

        if (!params.name) {
            log.error "--name is required"
            exit 1
        }
        if (params.fastq2) {
            Channel
                .fromList([tuple(params.name, file(params.fastq1), file(params.fastq2))])
                .set { input_fastqs }
        }
        else {
            Channel
                .fromList([tuple(params.name, file(params.fastq1))])
                .set { input_fastqs }
        }
    }
}
else if (params.input_fastas_list || params.fasta) {
    if (params.input_fastas_list) {
        Channel
            .fromPath(params.input_fastas_list)
            .splitCsv(header: ['name', 'fasta'], sep: "\t")
            .map{ row-> tuple(row.name, file(row.fasta)) }
            .set { input_fastas }
    }
    else {

        if (!params.name) {
            log.error "--name is required"
            exit 1
        }
        Channel
            .fromList([tuple(params.name, file(params.fasta))])
            .set { input_fastas }
    }
}
else {
    log.error "missing some input data"
    exit 1
}
if (params.skip_bcftools && params.skip_gatk && params.skip_ivar && params.skip_lofreq) {
    log.error "enable at least one variant caller"
    exit 1
}


workflow {

    if (input_fastqs) {
        if (params.library == "paired") {
            readTrimmingPairedEnd(input_fastqs)
            alignmentPairedEnd(readTrimmingPairedEnd.out[0], params.reference)
            bam_files = alignmentPairedEnd.out
        }
        else {
            readTrimmingSingleEnd(input_fastqs)
            alignmentSingleEnd(readTrimmingSingleEnd.out[0], params.reference)
            bam_files = alignmentSingleEnd.out
        }
        bamPreprocessing(bam_files, params.reference)

        // variant calling
        vcfs_to_normalize = null
        if (!params.skip_bcftools) {
            variantCallingBcfTools(bamPreprocessing.out[0], params.reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                variantCallingBcfTools.out : vcfs_to_normalize.concat(variantCallingBcfTools.out)
        }
        if (!params.skip_lofreq) {
            variantCallingLofreq(bamPreprocessing.out[0], params.reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                variantCallingLofreq.out : vcfs_to_normalize.concat(variantCallingLofreq.out)
        }
        if (!params.skip_gatk) {
            variantCallingGatk(bamPreprocessing.out[0], params.reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                variantCallingGatk.out : vcfs_to_normalize.concat(variantCallingGatk.out)
        }
        if (!params.skip_ivar) {
            variantCallingIvar(bamPreprocessing.out[0], params.reference, gff)
            // TODO: transform iVar to a VCF and follow normalization...
        }
    }
    else if (input_fastas) {
        // assembly variant calling
        assemblyVariantCaller(input_fastas, params.reference)
        vcfs_to_normalize = assemblyVariantCaller.out
    }

    variantNormalization(vcfs_to_normalize, params.reference)

    if (params.skip_sarscov2_annotations) {
        variantAnnotation(variantNormalization.out)
    }
    else {
        variantSarsCov2Annotation(variantNormalization.out)
    }
}
