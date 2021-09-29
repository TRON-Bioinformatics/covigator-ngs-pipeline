#!/usr/bin/env nextflow

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
else {
    reference = file(params.reference)
}

if (!params.gff && !params.fasta) {
    log.error "--gff is required"
    exit 1
}
else {
    gff = file(params.gff)
}

if (!params.name) {
    log.error "--name is required"
    exit 1
}

if (!params.fastq1 && !params.fasta) {
    log.error "either --fastq1 or --fasta are required"
    exit 1
}
else if (params.fastq1 && params.fasta) {
    log.error "provide only --fastq1 or --fasta"
    exit 1
}

if (params.skip_bcftools && params.skip_gatk && params.skip_ivar && params.skip_lofreq) {
    log.error "enable at least one variant caller"
    exit 1
}

if (!params.skip_sarscov2_annotations) {
    conservation_sarscov2 = file(params.conservation_sarscov2)
    conservation_sarscov2_header = file(params.conservation_sarscov2_header)
    conservation_sarbecovirus = file(params.conservation_sarbecovirus)
    conservation_sarbecovirus_header = file(params.conservation_sarbecovirus_header)
    conservation_vertebrate = file(params.conservation_vertebrate)
    conservation_vertebrate_header = file(params.conservation_vertebrate_header)
    pfam_names = file(params.pfam_names)
    pfam_descriptions = file(params.pfam_descriptions)
    pfam_names_header = file(params.pfam_names_header)
    pfam_descriptions_header = file(params.pfam_descriptions_header)
}


library = "paired"
if (!params.fastq2) {
    library = "single"
}

if (params.fastq1) {
    if (library == "paired") {

        process readTrimmingPairedEnd {
            cpus params.cpus
            memory params.memory
            tag params.name
            publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"

            input:
                val name from params.name
                file fastq1 from file(params.fastq1)
                file fastq2 from file(params.fastq2)

            output:
                set name, file("${fastq1.baseName}.trimmed.fq.gz"), file("${fastq2.baseName}.trimmed.fq.gz") into trimmed_fastqs
                file("${name}.fastp_stats.json")
                file("${name}.fastp_stats.html")

            """
            # --input_files needs to be forced, otherwise it is inherited from profile in tests
            fastp \
            --in1 ${fastq1} \
            --in2 ${fastq2} \
            --out1 ${fastq1.baseName}.trimmed.fq.gz \
            --out2 ${fastq2.baseName}.trimmed.fq.gz \
            --json ${name}.fastp_stats.json \
            --html ${name}.fastp_stats.html
            """
        }

        process alignmentPairedEnd {
            cpus params.cpus
            memory params.memory
            tag params.name

            input:
                set name, file(fastq1), file(fastq2) from trimmed_fastqs

            output:
                set name, file("${name}.bam") into bam_files

            """
            bwa mem -t ${task.cpus} ${reference} ${fastq1} ${fastq2} | \
            samtools view -uS - | \
            samtools sort - > ${name}.bam
            """
        }
    }
    else {

        process readTrimmingSingleEnd {
            cpus params.cpus
            memory params.memory
            tag params.name
            publishDir "${params.output}", mode: "copy", pattern: "*fastp_stats*"

            input:
                val name from params.name
                file fastq1 from file(params.fastq1)

            output:
                set name, file("${fastq1.baseName}.trimmed.fq.gz") into trimmed_fastqs
                file("${name}.fastp_stats.json")
                file("${name}.fastp_stats.html")

            """
            # --input_files needs to be forced, otherwise it is inherited from profile in tests
            fastp \
            --in1 ${fastq1} \
            --out1 ${fastq1.baseName}.trimmed.fq.gz \
            --json ${name}.fastp_stats.json \
            --html ${name}.fastp_stats.html
            """
        }

        process alignmentSingleEnd {
            cpus params.cpus
            memory params.memory
            tag params.name

            input:
                set name, file(fastq1) from trimmed_fastqs

            output:
                set name, file("${name}.bam") into bam_files

            """
            bwa mem -t ${task.cpus} ${reference} ${fastq1} | \
            samtools view -uS - | \
            samtools sort - > ${name}.bam
            """
        }
    }

    process bamPreprocessing {
        cpus params.cpus
        memory params.memory
        tag params.name
        if (params.keep_intermediate) {
            publishDir "${params.output}", mode: "copy"
        }
        publishDir "${params.output}", mode: "copy", pattern: "${name}.deduplication_metrics.txt"
        publishDir "${params.output}", mode: "copy", pattern: "${name}.coverage.tsv"
        publishDir "${params.output}", mode: "copy", pattern: "${name}.depth.tsv"

        input:
            set name, file(bam) from bam_files

        output:
            set name, file("${name}.preprocessed.bam"), file("${name}.preprocessed.bai") into preprocessed_bams,
                preprocessed_bams2, preprocessed_bams3, preprocessed_bams4
            file "${name}.deduplication_metrics.txt"
            file "${name}.coverage.tsv"
            file "${name}.depth.tsv"


        """
        gatk CleanSam \
        --java-options '-Xmx${params.memory} -Djava.io.tmpdir=tmp' \
        --INPUT ${bam} \
        --OUTPUT /dev/stdout | \
        gatk AddOrReplaceReadGroups \
        --java-options '-Xmx${params.memory} -Djava.io.tmpdir=tmp' \
        --VALIDATION_STRINGENCY SILENT \
        --INPUT /dev/stdin \
        --OUTPUT ${bam.baseName}.prepared.bam \
        --REFERENCE_SEQUENCE ${reference} \
        --RGPU 1 \
        --RGID 1 \
        --RGSM ${name} \
        --RGLB 1 \
        --RGPL ILLUMINA \
        --SORT_ORDER queryname

        gatk MarkDuplicates \
        --java-options '-Xmx${params.memory}  -Djava.io.tmpdir=tmp' \
        --INPUT ${bam.baseName}.prepared.bam \
        --METRICS_FILE ${name}.deduplication_metrics.txt \
        --OUTPUT ${bam.baseName}.dedup.bam \
        --REMOVE_DUPLICATES true

        gatk SortSam \
        --java-options '-Xmx${params.memory}  -Djava.io.tmpdir=tmp' \
        --INPUT ${bam.baseName}.dedup.bam \
        --OUTPUT ${bam.baseName}.dedup.sorted.bam \
        --SORT_ORDER coordinate

        gatk BuildBamIndex --INPUT ${bam.baseName}.dedup.sorted.bam

        gatk3 -Xmx${params.memory} -Djava.io.tmpdir=tmp -T RealignerTargetCreator \
	    --input_file ${bam.baseName}.dedup.sorted.bam \
	    --out ${bam.baseName}.RA.intervals \
	    --reference_sequence ${reference}

	    gatk3 -Xmx${params.memory} -Djava.io.tmpdir=tmp -T IndelRealigner \
	    --input_file ${bam.baseName}.dedup.sorted.bam \
	    --out ${name}.preprocessed.bam \
	    --reference_sequence ${reference} \
	    --targetIntervals ${bam.baseName}.RA.intervals \
	    --consensusDeterminationModel USE_SW \
	    --LODThresholdForCleaning 0.4 \
	    --maxReadsInMemory 600000

	    samtools coverage ${name}.preprocessed.bam > ${name}.coverage.tsv

	    samtools depth -s -d 0 -H ${name}.preprocessed.bam > ${name}.depth.tsv
        """
    }

    vcfs_to_normalize = null

    if (!params.skip_bcftools) {
        process variantCallingBcfTools {
            cpus params.cpus
            memory params.memory
            tag params.name
            if (params.keep_intermediate) {
                publishDir "${params.output}", mode: "copy"
            }

            input:
                set name, file(bam), file(bai) from preprocessed_bams

            output:
                set name, file("${name}.bcftools.bcf") into bcftools_vcfs

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
             --ploidy 1 | \
            bcftools filter \
            --exclude 'INFO/IMF < ${params.low_frequency_variant_threshold}' \
            --soft-filter LOW_FREQUENCY - | \
            bcftools filter \
            --exclude 'INFO/IMF >= ${params.low_frequency_variant_threshold} && INFO/IMF < ${params.subclonal_variant_threshold}' \
            --soft-filter SUBCLONAL \
             --output-type b - > ${name}.bcftools.bcf
            """
        }
        vcfs_to_normalize = vcfs_to_normalize == null? bcftools_vcfs : vcfs_to_normalize.concat(bcftools_vcfs)
    }

    if (!params.skip_lofreq) {
        process variantCallingLofreq {
            cpus params.cpus
            memory params.memory
            tag params.name
            if (params.keep_intermediate) {
                publishDir "${params.output}", mode: "copy"
            }

            input:
                set name, file(bam), file(bai) from preprocessed_bams2

            output:
                set name, file("${name}.lofreq.vcf") into lofreq_vcfs

            """
            lofreq call \
            --min-bq ${params.min_base_quality} \
            --min-alt-bq ${params.min_base_quality} \
            --min-mq ${params.min_mapping_quality} \
            --ref ${reference} \
            --call-indels \
            <( lofreq indelqual --dindel --ref ${reference} ${bam} ) | \
            bgzip -c > ${name}.lofreq.vcf.gz

            tabix -p vcf ${name}.lofreq.vcf.gz

            # annotates low frequency and subclonal variants
            bcftools view -Ob ${name}.lofreq.vcf.gz | \
            bcftools filter \
            --exclude 'INFO/AF < ${params.low_frequency_variant_threshold}' \
            --soft-filter LOW_FREQUENCY - | \
            bcftools filter \
            --exclude 'INFO/AF >= ${params.low_frequency_variant_threshold} && INFO/AF < ${params.subclonal_variant_threshold}' \
            --soft-filter SUBCLONAL - > ${name}.lofreq.vcf
            """
        }
        vcfs_to_normalize = vcfs_to_normalize == null? lofreq_vcfs : vcfs_to_normalize.concat(lofreq_vcfs)
    }

    if (!params.skip_gatk) {
        process variantCallingGatk {
            cpus params.cpus
            memory params.memory
            tag params.name
            if (params.keep_intermediate) {
                publishDir "${params.output}", mode: "copy"
            }

            input:
                set name, file(bam), file(bai) from preprocessed_bams3

            output:
                set name, file("${name}.gatk.vcf") into gatk_vcfs

            """
            gatk HaplotypeCaller \
            --input $bam \
            --output ${name}.gatk.vcf \
            --reference ${reference} \
            --ploidy 1 \
            --min-base-quality-score ${params.min_base_quality} \
            --minimum-mapping-quality ${params.min_mapping_quality} \
            --annotation AlleleFraction
            """
        }
        vcfs_to_normalize = vcfs_to_normalize == null? gatk_vcfs : vcfs_to_normalize.concat(gatk_vcfs)
    }

    if (!params.skip_ivar) {
        process variantCallingIvar {
            cpus params.cpus
            memory params.memory
            tag params.name
            publishDir "${params.output}", mode: "copy"

            input:
                set name, file(bam), file(bai) from preprocessed_bams4

            output:
                file("${name}.ivar.tsv")

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
    }
}
else if (params.fasta) {

    process assemblyVariantCaller {
        cpus params.cpus
        memory params.memory
        tag params.name
        if (params.keep_intermediate) {
            publishDir "${params.output}", mode: "copy"
        }

        input:
            val name from params.name
            file fasta from file(params.fasta)

        output:
            set name, file("${name}.assembly.vcf") into vcfs_to_normalize

        """
        assembly_variant_caller.py \
        --fasta ${fasta} \
        --reference ${reference} \
        --output-vcf ${name}.assembly.vcf \
        --match-score $params.match_score \
        --mismatch-score $params.mismatch_score \
        --open-gap-score $params.open_gap_score \
        --extend-gap-score $params.extend_gap_score \
        --chromosome $params.chromosome
        """
    }
}

process variantNormalization {
    cpus params.cpus
    memory params.memory
    tag params.name
    if (params.keep_intermediate) {
        publishDir "${params.output}", mode: "copy"
    }

    input:
        set name, file(vcf) from vcfs_to_normalize

    output:
      set name, file("${vcf.baseName}.normalized.vcf") into normalized_vcf_files

    script:
    """
    # initial sort of the VCF
    bcftools sort ${vcf} | \

    # checks reference genome, decompose multiallelics, trim and left align indels
    bcftools norm --multiallelics -any --check-ref e --fasta-ref ${params.reference} \
    --old-rec-tag OLD_CLUMPED - | \

    # decompose complex variants
    vt decompose_blocksub -a -p - | \

    # remove duplicates after normalisation
    bcftools norm --rm-dup exact -o ${vcf.baseName}.normalized.vcf -
    """
}

if (params.skip_sarscov2_annotations) {
    process variantAnnotation {
        cpus params.cpus
        memory params.memory
        tag params.name
        publishDir "${params.output}", mode: "copy"

        input:
            set name, file(vcf) from normalized_vcf_files

        output:
            file("${vcf.baseName}.annotated.vcf.gz")
            file("${vcf.baseName}.annotated.vcf.gz.tbi")

        """
        # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
        cp ${params.snpeff_config} .

        snpEff eff -dataDir ${params.snpeff_data} \
        -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
        ${params.snpeff_organism}  ${vcf} | \
        bgzip -c > ${vcf.baseName}.annotated.vcf.gz

        tabix -p vcf ${vcf.baseName}.annotated.vcf.gz
        """
    }
} else {
    process variantSarsCov2Annotation {
        cpus params.cpus
        memory params.memory
        tag params.name
        publishDir "${params.output}", mode: "copy"

        input:
            set name, file(vcf) from normalized_vcf_files

        output:
            file("${vcf.baseName}.annotated.vcf.gz")
            file("${vcf.baseName}.annotated.vcf.gz.tbi")

        """
        # for some reason the snpEff.config file needs to be in the folder where snpeff runs...
        cp ${params.snpeff_config} .

        snpEff eff -dataDir ${params.snpeff_data} \
        -noStats -no-downstream -no-upstream -no-intergenic -no-intron -onlyProtein -hgvs1LetterAa -noShiftHgvs \
        Sars_cov_2.ASM985889v3.101  ${vcf} | \
        bgzip -c | \
        bcftools annotate \
        --annotations ${conservation_sarscov2} \
        --header-lines ${conservation_sarscov2_header} \
        -c CHROM,FROM,TO,CONS_HMM_SARS_COV_2 \
        --output-type z - | \
        bcftools annotate \
        --annotations ${conservation_sarbecovirus} \
        --header-lines ${conservation_sarbecovirus_header} \
        -c CHROM,FROM,TO,CONS_HMM_SARBECOVIRUS \
        --output-type z - | \
        bcftools annotate \
        --annotations ${conservation_vertebrate} \
        --header-lines ${conservation_vertebrate_header} \
        -c CHROM,FROM,TO,CONS_HMM_VERTEBRATE_COV \
        --output-type z - | \
        bcftools annotate \
        --annotations ${pfam_names} \
        --header-lines ${pfam_names_header} \
        -c CHROM,FROM,TO,PFAM_NAME \
        --output-type z - | \
        bcftools annotate \
        --annotations ${pfam_descriptions} \
        --header-lines ${pfam_descriptions_header} \
        -c CHROM,FROM,TO,PFAM_DESCRIPTION \
        --output-type z - > ${vcf.baseName}.annotated.vcf.gz

        # TODO: include this step for GISAID data
        #bcftools annotate \
        #--annotations ${params.problematic_sites} \
        #--columns FILTER \
        #--output-type b - > ${vcf.baseName}.annotated.vcf.gz

        tabix -p vcf ${vcf.baseName}.annotated.vcf.gz
        """
    }
}
