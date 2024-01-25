#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { READ_TRIMMING_PAIRED_END; READ_TRIMMING_SINGLE_END } from './modules/01_fastp'
include { ALIGNMENT_PAIRED_END; ALIGNMENT_SINGLE_END } from './modules/02_bwa'
include { BAM_PREPROCESSING; COVERAGE_ANALYSIS; PRIMER_TRIMMING_IVAR } from './modules/03_bam_preprocessing'
include { VARIANT_CALLING_BCFTOOLS; VARIANT_CALLING_LOFREQ ; VARIANT_CALLING_GATK ;
            VARIANT_CALLING_IVAR ; VARIANT_CALLING_ASSEMBLY; IVAR2VCF } from './modules/04_variant_calling'
include { VARIANT_NORMALIZATION ; PHASING } from './modules/05_variant_normalization'
include { VARIANT_ANNOTATION; VARIANT_SARSCOV2_ANNOTATION;
            VARIANT_VAF_ANNOTATION; VAFATOR } from './modules/06_variant_annotation'
include { PANGOLIN_LINEAGE; VCF2FASTA } from './modules/07_lineage_annotation'
include { BGZIP } from './modules/08_compress_vcf'

params.help= false

params.fastq1 = false
params.fastq2 = false
params.fasta = false
params.vcf = false
params.bam = false
params.name = false

params.skip_lofreq = false
params.skip_ivar = false
params.skip_bcftools = false
params.skip_gatk = false
params.skip_pangolin = false
params.skip_normalization = false

// references
params.reference = false
params.gff = false
params.snpeff_data = false
params.snpeff_config = false
params.snpeff_organism = false
params.primers = false

params.output = "."
params.min_mapping_quality = 20
params.min_base_quality = 20
params.vafator_min_mapping_quality = 0
params.vafator_min_base_quality = 0
params.low_frequency_variant_threshold = 0.02
params.subclonal_variant_threshold = 0.5
params.lq_clonal_variant_threshold = 0.8
params.memory = "3g"
params.cpus = 1
params.keep_intermediate = false
params.match_score = 2
params.mismatch_score = -1
params.open_gap_score = -3
params.extend_gap_score = -0.1
params.skip_sarscov2_annotations = false
params.library = false
params.input_fastqs_list = false
params.input_fastas_list = false
params.input_vcfs_list = false
params.input_bams_list = false

if (params.help) {
    log.info params.help_message
    exit 0
}
if (params.output == false) {
    log.error "--output is required"
    exit 1
}
if (params.fastq1 != false && params.fasta != false) {
    log.error "provide only --fastq1 or --fasta"
    exit 1
}
if (params.input_fastqs_list != false && params.input_fastas_list != false) {
    log.error "provide only --input_fastqs_list or --input_fastas_list"
    exit 1
}

if (params.reference == false) {
    log.info "Using default SARS-CoV-2 reference genome"
    reference = params.sarscov2_reference   // do not put into a file as we need the indices
    gff = file(params.sarscov2_gff)
    snpeff_data = params.sarscov2_snpeff_data
    snpeff_config = params.sarscov2_snpeff_config
    snpeff_organism = params.sarscov2_snpeff_organism
    skip_sarscov2_annotations = params.skip_sarscov2_annotations
}
else {
    log.info "Using custom reference genome: ${params.reference}"
    reference = params.reference    // do not put into a file as we need the indices
    gff = params.gff ? file(params.gff) : false
    snpeff_data = params.snpeff_data
    snpeff_config = params.snpeff_config
    snpeff_organism = params.snpeff_organism
    skip_sarscov2_annotations = true
}

primers = params.primers ? file(params.primers) : false

skip_snpeff = false
if (! snpeff_data || ! snpeff_config || ! snpeff_organism) {
    log.info "Skipping SnpEff annotation as either --snpeff_data, --snpeff_config or --snpeff_organism was not provided"
    skip_snpeff = true
}

input_fastqs = false
input_fastas = false
input_vcfs = false
preprocessed_bams = false
library = params.library
if (params.input_fastqs_list != false || params.fastq1 != false) {

    // if independent FASTQ files are provided the value of library is overridden
    if (params.fastq1 != false && params.fastq2 == false) {
        library = "single"
    }
    else if (params.fastq1 != false && params.fastq2 != false) {
        library = "paired"
    }
    else if (params.input_fastqs_list && library == false) {
        log.error "--library paired|single is required when --input_fastqs_list is provided"
        exit 1
    }

    if (params.input_fastqs_list) {
        if (library == "paired") {
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

        if (params.name == false) {
            log.error "--name is required"
            exit 1
        }
        if (params.fastq2 != false) {
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
            .map{ row-> tuple(row.name, "assembly", file(row.fasta)) }
            .set { input_fastas }
    }
    else {

        if (params.name == false) {
            log.error "--name is required"
            exit 1
        }
        Channel
            .fromList([tuple(params.name, "assembly", file(params.fasta))])
            .set { input_fastas }
    }
}
else if (params.input_vcfs_list != false || params.vcf != false) {
    if (params.input_vcfs_list) {
        Channel
            .fromPath(params.input_vcfs_list)
            .splitCsv(header: ['name', 'vcf'], sep: "\t")
            .map{ row-> tuple(row.name, "input", file(row.vcf)) }
            .set { input_vcfs }
    }
    else {

        if (params.name == false) {
            log.error "--name is required"
            exit 1
        }
        Channel
            .fromList([tuple(params.name, "input", file(params.vcf))])
            .set { input_vcfs }
    }
}
else {
    log.error "missing some input data"
    exit 1
}

if (params.input_bams_list || params.bam) {
    if (params.input_bams_list) {
        Channel
            .fromPath(params.input_bams_list)
            .splitCsv(header: ['name', 'bam', 'bai'], sep: "\t")
            .map{ row-> tuple(row.name, file(row.bam), file(row.bai)) }
            .set { preprocessed_bams }
    }
    else {
        if (params.name == false) {
            log.error "--name is required"
            exit 1
        }
        if (params.bai == false) {
            log.error "--bai is required"
            exit 1
        }
        Channel
            .fromList([tuple(params.name, file(params.bam), file(params.bai))])
            .set { preprocessed_bams }
    }
}


if (params.skip_bcftools && params.skip_gatk && params.skip_ivar && params.skip_lofreq) {
    log.error "enable at least one variant caller"
    exit 1
}


workflow {
    if (input_fastqs) {
        if (library == "paired") {
            READ_TRIMMING_PAIRED_END(input_fastqs)
            ALIGNMENT_PAIRED_END(READ_TRIMMING_PAIRED_END.out[0], reference)
            bam_files = ALIGNMENT_PAIRED_END.out
        }
        else {
            READ_TRIMMING_SINGLE_END(input_fastqs)
            ALIGNMENT_SINGLE_END(READ_TRIMMING_SINGLE_END.out[0], reference)
            bam_files = ALIGNMENT_SINGLE_END.out
        }
        BAM_PREPROCESSING(bam_files, reference)
        preprocessed_bams = BAM_PREPROCESSING.out.preprocessed_bams

        if (primers) {
            PRIMER_TRIMMING_IVAR(preprocessed_bams, primers)
            preprocessed_bams = PRIMER_TRIMMING_IVAR.out.trimmed_bam
        }
        COVERAGE_ANALYSIS(preprocessed_bams)

        // variant calling
        vcfs_to_normalize = null
        if (!params.skip_bcftools) {
            VARIANT_CALLING_BCFTOOLS(preprocessed_bams, reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                VARIANT_CALLING_BCFTOOLS.out : vcfs_to_normalize.concat(VARIANT_CALLING_BCFTOOLS.out)
        }
        if (!params.skip_lofreq) {
            VARIANT_CALLING_LOFREQ(preprocessed_bams, reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                VARIANT_CALLING_LOFREQ.out : vcfs_to_normalize.concat(VARIANT_CALLING_LOFREQ.out)
        }
        if (!params.skip_gatk) {
            VARIANT_CALLING_GATK(preprocessed_bams, reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                VARIANT_CALLING_GATK.out : vcfs_to_normalize.concat(VARIANT_CALLING_GATK.out)
        }
        if (!params.skip_ivar && gff) {
            VARIANT_CALLING_IVAR(preprocessed_bams, reference, gff)
            IVAR2VCF(VARIANT_CALLING_IVAR.out, reference)
            vcfs_to_normalize = vcfs_to_normalize == null?
                IVAR2VCF.out : vcfs_to_normalize.concat(IVAR2VCF.out)
        }
    }
    else if (input_fastas) {
        if (!params.skip_pangolin) {
            // pangolin from fasta
            PANGOLIN_LINEAGE(input_fastas)
        }

        // assembly variant calling
        VARIANT_CALLING_ASSEMBLY(input_fastas, reference)
        vcfs_to_normalize = VARIANT_CALLING_ASSEMBLY.out
    }
    else if (input_vcfs) {
        vcfs_to_normalize = input_vcfs
    }

    if (! params.skip_normalization) {
        VARIANT_NORMALIZATION(vcfs_to_normalize, reference)
        normalized_vcfs = VARIANT_NORMALIZATION.out
    }
    else {
        normalized_vcfs = vcfs_to_normalize
    }

    if (input_fastqs || input_vcfs) {
        // pangolin from VCF on the normalized VCFs
        if (!params.skip_pangolin) {
            VCF2FASTA(normalized_vcfs, reference)
            PANGOLIN_LINEAGE(VCF2FASTA.out)
        }
    }

    if (! skip_sarscov2_annotations) {
        // only optionally add SARS-CoV-2 specific annotations
        VARIANT_SARSCOV2_ANNOTATION(normalized_vcfs)
        normalized_vcfs = VARIANT_SARSCOV2_ANNOTATION.out.annotated_vcfs
    }

    if (preprocessed_bams) {
        // we can only add technical annotations when we have the reads
        VAFATOR(normalized_vcfs.combine(preprocessed_bams, by: 0))
        VARIANT_VAF_ANNOTATION(VAFATOR.out.annotated_vcf)
        normalized_vcfs = VARIANT_VAF_ANNOTATION.out.vaf_annotated
    }

    // NOTE: phasing has to happen before SnpEff annotation for MNVs to be annotated correctly
    if (gff) {
        PHASING(normalized_vcfs, reference, gff)
        normalized_vcfs = PHASING.out
    }

    if (! skip_snpeff) {
        // only when configured we run SnpEff
        VARIANT_ANNOTATION(normalized_vcfs, snpeff_data, snpeff_config, snpeff_organism)
        normalized_vcfs = VARIANT_ANNOTATION.out.annotated_vcfs
    }
    else {
        BGZIP(normalized_vcfs)
    }
}
