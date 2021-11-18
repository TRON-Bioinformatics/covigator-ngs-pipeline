![CoVigator logo](images/CoVigator_logo_txt_nobg.png "CoVigator logo")

# Covigator NGS pipeline: full variant detection pipeline for Sars-CoV-2

[![DOI](https://zenodo.org/badge/374669617.svg)](https://zenodo.org/badge/latestdoi/374669617)
[![Run tests](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline/actions/workflows/automated_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/covigator-ngs-pipeline/actions/workflows/automated_tests.yml)
[![Powered by NumFOCUS](https://img.shields.io/badge/powered%20by-Nextflow-orange.svg?style=flat&colorA=E1523D&colorB=007D8A)](https://www.nextflow.io/)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)



The Covigator pipeline processes SARS-CoV-2 FASTQ or FASTA files into annotated and normalized analysis ready VCF files. 
The pipeline is implemented in the Nextflow framework (Di Tommaso, 2017).

## Possible inputs

- FASTQ files, either one FASTQ file for single end or two FASTQ files for paired end

or

- FASTA file with an assembly (ie: a single DNA sequence)

## Outputs

- Multiple VCF files from different variant callers when FASTQ files are provided

or

- A single VCF file from the global alignment of the assembly against the reference genome


## Pipeline details

When FASTQ files are provided the pipeline includes the following steps:
- **Trimming**. `fastp` is used to trim reads with default values. This step also includes QC filtering.
- **Alignment**. `BWA mem` is used for the alignment of single or paired end samples.
- **BAM preprocessing**. BAM files are prepared and duplicate reads are marked using GATK and Picard tools.
- **Coverage analysis**. `samtools coverage` and `samtools depth` are used to compute the horizontal and vertical 
  coverage respectively.
- **Variant calling**. Four different variant callers are employed: BCFtools, LoFreq, iVar and GATK. 
  Subsequent processing of resulting VCF files is independent for each caller, except for iVar which does not produce a VCF file but a custom TSV file.
- **Variant normalization**. `bcftools norm` and `vt` tools are employed to left align indels, trim variant calls and remove variant duplicates.
- **Variant annotation**. `SnpEff` is employed to annotate the variant consequences of variants, 
  `bcftools annotate` is employed to add additional annotations.

Both single end and paired end FASTQ files are supported.

When a FASTA file is provided with a single assembly sequence the pipeline includes the following steps:
- **Variant calling**. A Smith-Waterman global alignment is performed against the reference sequence to call SNVs and 
  indels. Indels longer than 50 bp and at the beginning or end of the assembly sequence are excluded. Any mutation where
  either reference or assembly contain a N is excluded.
- **Variant normalization**. Same as described above.
- **Variant annotation**. Same as described above.

The FASTA file is expected to contain a single assembly sequence. 
Bear in mind that only clonal variants can be called on the assembly.

The alignment, BAM preprocessing and variant normalization pipelines are based on the implementations in additional 
Nextflow pipelines within the TronFlow initiative. 
The full details are available in their respective repositories:
- https://github.com/TRON-Bioinformatics/tronflow-bwa (https://doi.org/10.5281/zenodo.4722852)
- https://github.com/TRON-Bioinformatics/tronflow-bam-preprocessing (https://doi.org/10.5281/zenodo.4810918)
- https://github.com/TRON-Bioinformatics/tronflow-variant-normalization (https://doi.org/10.5281/zenodo.4875095)

The default SARS-CoV-2 reference files correspond to Sars_cov_2.ASM985889v3 and were downloaded from Ensembl servers.
These references can be customised to use a different SARS-CoV-2 reference or to analyse a different virus.
Two files need to be provided: a sequence file in FASTA format and a gene annotation file in GFFv3 format. 
Additionally, the FASTA needs bwa indexes and .fai index.
These indexes can be generated with the following two commands:
```
bwa index reference.fasta
samtools faidx reference.fasta
```

The LoFreq variants are annotated on the `FILTER` column using the reported variant allele frequency 
(VAF) into `LOW_FREQUENCY`, `SUBCLONAL` and finally `PASS` variants correspond to clonal variants. By default, 
variants with a VAF < 20 % are considered `LOW_FREQUENCY` and variants with a VAF >= 20 % and < 80 % are considered 
`SUBCLONAL`. This thresholds can be changed with the parameters `--low_frequency_variant_threshold` and
`--subclonal_variant_threshold`. Indels called by BCFtools are also annotated by VAF, but not SNVs.

All variant calls are additionally annotated with:
- ConsHMM conservation scores as reported in (Kwon, 2021)
- Pfam domains as reported in Ensemble annotations.

A variant in the output VCF will look as follows:
```
MN908947.3      21680   .       G       A       250     LOW_FREQUENCY   DP=1252;AF=0.015176;SB=4;DP4=551,679,11,8;ANN=A|missense_variant|MODERATE|S|gene-GU280_gp02|transcript|TRANSCRIPT_gene-GU280_gp02|protein_coding|1/1|c.118G>A|p.D40N|118/3822|118/3822|40/1273||;CONS_HMM_SARS_COV_2=0.57215;CONS_HMM_SARBECOVIRUS=0.57215;CONS_HMM_VERTEBRATE_COV=0;PFAM_NAME=bCoV_S1_N;PFAM_DESCRIPTION=Betacoronavirus-like spike glycoprotein S1, N-terminal
```

Where:
- `INFO/DP` is the number of reads overlapping this position
- `INFO/AF` is the VAF as reported by LoFreq (only avalable in LoFreq calls)
- `INFO/ANN` are the SnpEff consequence annotations
- `INFO/CONS_HMM_SARS_COV_2` is the ConsHMM conservation score in SARS-CoV-2
- `INFO/CONS_HMM_SARBECOVIRUS` is the ConsHMM conservation score among Sarbecovirus
- `INFO/CONS_HMM_VERTEBRATE_COV` is the ConsHMM conservation score among vertebrate Corona virus
- `INFO/PFAM_NAME` is the Interpro name for the overlapping Pfam domains
- `INFO/PFAM_DESCRIPTION` is the Interpro description for the overlapping Pfam domains


## Requirements

- Nextflow >= 19.10.0
- Java >= 8
- Conda >=4.9

## How to run it

For paired end reads:
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] --fastq1 <FASTQ_FILE> --fastq2 <FASTQ_FILE> --name example_run --output <OUTPUT_FOLDER> [--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] [--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

For single end reads:
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] --fastq1 <FASTQ_FILE> --name example_run --output <OUTPUT_FOLDER> [--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] [--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

For assembly:
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] --fasta <FASTA_FILE> --name example_run --output <OUTPUT_FOLDER> [--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] [--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```

For batch processing of reads use `--input_fastqs_list` and `--name`.
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] --input_fastqs_list <TSV_FILE> --library <paired|single> --output <OUTPUT_FOLDER> [--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] [--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```
where the TSV file contains two or three columns tab-separated columns without header. Columns: sample name, path to FASTQ 1 and optionally path to FASTQ 2. 

For batch processing of assemblies use `--input_fastas_list`.
```
nextflow run tron-bioinformatics/covigator-ngs-pipeline [-profile conda] --input_fastas_list <TSV_FILE> --library <paired|single> --output <OUTPUT_FOLDER> [--reference <path_to_reference>/Sars_cov_2.ASM985889v3.fa] [--gff <path_to_reference>/Sars_cov_2.ASM985889v3.gff3]
```
where the TSV file contains two columns tab-separated columns without header. Columns: sample name and path to FASTA.

All options available using `--help`:
```
$ nextflow run tron-bioinformatics/covigator-ngs-pipeline -profile conda --help

Usage:
    nextflow run tron-bioinformatics/covigator-ngs-pipeline -profile conda --help

Input:
    * --fastq1: the first input FASTQ file (not compatible with --fasta)
    * --fasta: the FASTA file containing the assembly sequence (not compatible with --fastq1)
    * --name: the sample name, output files will be named after this name
    * --reference: the reference genome FASTA file, *.fai, *.dict and bwa indexes are required.
    * --gff: the GFFv3 gene annotations file (only required with --fastq1)
    * --output: the folder where to publish output
    * --input_fastqs_list: alternative to --name and --fastq1 for batch processing
    * --library: required only when using --input_fastqs
    * --input_fastas_list: alternative to --name and --fasta for batch processing

Optional input:
    * --fastq2: the second input FASTQ file
    * --min_base_quality: minimum base call quality to take a base into account (default: 20)
    * --min_mapping_quality: minimum mapping quality to take a read into account (default: 20)
    * --low_frequency_variant_threshold: VAF threshold to mark a variant as low frequency (default: 0.2)
    * --subclonal_variant_threshold: VAF superior threshold to mark a variant as subclonal (default: 0.8)
    * --memory: the ammount of memory used by each job (default: 3g)
    * --cpus: the number of CPUs used by each job (default: 1)
    * --initialize: initialize the conda environment
    * --skip_lofreq: skips calling variants with LoFreq
    * --skip_gatk: skips calling variants with GATK
    * --skip_bcftools: skips calling variants with BCFTools
    * --skip_ivar: skips calling variants with iVar
    * --match_score: global alignment match score, only applicable for assemblies (default: 2)
    * --mismatch_score: global alignment mismatch score, only applicable for assemblies (default: -1)
    * --open_gap_score: global alignment open gap score, only applicable for assemblies (default: -3)
    * --extend_gap_score: global alignment extend gap score, only applicable for assemblies (default: -0.1)
    * --chromosome: chromosome for variant calls, only applicable for assemblies (default: "MN908947.3")
    * --skip_sarscov2_annotations: skip some of the SARS-CoV-2 specific annotations (default: false)
    * --snpeff_data: path to the SnpEff data folder, it will be useful to use the pipeline on other virus than SARS-CoV-2
    * --snpeff_config: path to the SnpEff config file, it will be useful to use the pipeline on other virus than SARS-CoV-2
    * --snpeff_organism: organism to annotate with SnpEff, it will be useful to use the pipeline on other virus than SARS-CoV-2

Output:
    * Output a normalized, phased and annotated VCF file for each of BCFtools, GATK and LoFreq when FASTQ files are
    provided or a single VCF obtained from a global alignment when a FASTA file is provided
    * Output a TSV file output from iVar
```

### Initializing the conda environments

If you are planning to use it concurrently on multiple samples with conda, first initialize the conda environment, 
otherwise concurrent sample executions will clash trying to create the same conda environment. Run:
```
nextflow main.nf -profile conda --initialize
```

This will create the necessary conda environment under `work/conda`, subsequent concurrent executions will use the same
conda environment.


## References

- Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. Unified Representation of Genetic Variants. Bioinformatics (2015) 31(13): 2202-2204](http://bioinformatics.oxfordjournals.org/content/31/13/2202) and uses bcftools [Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (Oxford, England), 27(21), 2987–2993. 10.1093/bioinformatics/btr509
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.
- Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M. (2013). From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. Curr Protoc Bioinformatics, 43:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.
- Martin, M., Patterson, M., Garg, S., O Fischer, S., Pisanti, N., Klau, G., Schöenhuth, A., & Marschall, T. (2016). WhatsHap: fast and accurate read-based phasing. BioRxiv, 085050. https://doi.org/10.1101/085050
- Danecek, P., & McCarthy, S. A. (2017). BCFtools/csq: haplotype-aware variant consequences. Bioinformatics, 33(13), 2037–2039. https://doi.org/10.1093/bioinformatics/btx100
- Wilm, A., Aw, P. P. K., Bertrand, D., Yeo, G. H. T., Ong, S. H., Wong, C. H., Khor, C. C., Petric, R., Hibberd, M. L., & Nagarajan, N. (2012). LoFreq: A sequence-quality aware, ultra-sensitive variant caller for uncovering cell-population heterogeneity from high-throughput sequencing datasets. Nucleic Acids Research, 40(22), 11189–11201. https://doi.org/10.1093/nar/gks918
- Grubaugh, N. D., Gangavarapu, K., Quick, J., Matteson, N. L., De Jesus, J. G., Main, B. J., Tan, A. L., Paul, L. M., Brackney, D. E., Grewal, S., Gurfield, N., Van Rompay, K. K. A., Isern, S., Michael, S. F., Coffey, L. L., Loman, N. J., & Andersen, K. G. (2019). An amplicon-based sequencing framework for accurately measuring intrahost virus diversity using PrimalSeq and iVar. Genome Biology, 20(1), 8. https://doi.org/10.1186/s13059-018-1618-7
- Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu; fastp: an ultra-fast all-in-one FASTQ preprocessor, Bioinformatics, Volume 34, Issue 17, 1 September 2018, Pages i884–i890, https://doi.org/10.1093/bioinformatics/bty560
- Kwon, S. Bin, & Ernst, J. (2021). Single-nucleotide conservation state annotation of the SARS-CoV-2 genome. Communications Biology, 4(1), 1–11. https://doi.org/10.1038/s42003-021-02231-w


## Resources

SARS-CoV-2 ASM985889v3 references were downloaded from Ensembl on 6th of October 2020:
- ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/dna/Sars_cov_2.ASM985889v3.dna.toplevel.fa.gz
- ftp://ftp.ensemblgenomes.org/pub/viruses/gff3/sars_cov_2/Sars_cov_2.ASM985889v3.101.gff3.gz

ConsHMM mutation depletion scores downloaded on 1st of July 2021:
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionConsHMM.bed
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionSarbecovirusConsHMM.bed
- https://github.com/ernstlab/ConsHMM_CoV/blob/master/wuhCor1.mutDepletionVertebrateCoVConsHMM.bed

Gene annotations including Pfam domains downloaded from Ensembl on 25th of February 2021 from:
- ftp://ftp.ensemblgenomes.org/pub/viruses/json/sars_cov_2/sars_cov_2.json