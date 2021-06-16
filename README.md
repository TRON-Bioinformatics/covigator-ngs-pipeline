# Covigator NGS pipeline

[![DOI](https://zenodo.org/badge/374669617.svg)](https://zenodo.org/badge/latestdoi/374669617)

The Covigator NGS pipeline process SARS-CoV-2 FASTQ files into analysis ready VCF files. The pipeline is implemented in the Nextflow framework (Di Tommaso, 2017).

The pipeline includes the following steps:
- **Trimming**. `fastp` is used to trim reads with default values.
- **Alignment**. `BWA mem` is used for the alignment of single or paired end samples.
- **BAM preprocessing**. BAM files are prepared and duplicate reads are marked using GATK and Picard tools.
- **Variant calling**. Four different variant callers are employed: BCFtools, LoFreq, iVar and GATK. 
  Subsequent processing of resulting VCF files is independent for each caller, except for iVar which does not produce a VCF file but a custom TSV file.
- **Variant normalization**. `bcftools norm` and `vt` tools are employed to left align indels, trim variant calls and remove variant duplicates.
- **Variant consequence annotation**. `SnpEff` is employed to annotate the variant consequences of variants.

The alignment, BAM preprocessing and variant normalization pipelines were implemented in additional Nextflow pipelines within the TronFlow initiative. 
The full details are available in their respective repositories:
- https://github.com/TRON-Bioinformatics/tronflow-bwa (https://doi.org/10.5281/zenodo.4722852)
- https://github.com/TRON-Bioinformatics/tronflow-bam-preprocessing (https://doi.org/10.5281/zenodo.4810918)
- https://github.com/TRON-Bioinformatics/tronflow-variant-normalization (https://doi.org/10.5281/zenodo.4875095)

The default SARS-CoV-2 reference files correspond to Sars_cov_2.ASM985889v3 and were downloaded from Ensembl servers.
These references can be customised to use a different SARS-CoV-2 reference or to analyse a different virus.
Two files need to be provided FASTA and GFFv3. Additionally, the FASTA needs several indexes: bwa indexes, .fai index and .dict index.
These indexes can be generated with the following commands:
```
bwa index reference.fasta
samtools faidx reference.fasta
gatk CreateSequenceDictionary -R reference.fasta
```

The LoFreq variants are annotated on the `FILTER` column using the reported variant allele frequency 
(VAF) into `LOW_FREQUENCY`, `SUBCLONAL` and finally `PASS` variants correspond to clonal variants. By default, 
variants with a VAF < 20 % are considered `LOW_FREQUENCY` and variants with a VAF >= 20 % and < 80 % are considered 
`SUBCLONAL`. This thresholds can be changed with the parameters `--low_frequency_variant_threshold` and
`--subclonal_variant_threshold`. Indels called by BCFtools are also annotated by VAF, but not SNVs.


## Requirements

- Nextflow >=19.10.0
- Java >= 8
- Conda >=4.9

## How to run it

```
$ nextflow run tron-bioinformatics/covigator-ngs-pipeline -profile conda --help

Usage:
    nextflow run tron-bioinformatics/covigator-ngs-pipeline -profile conda --help

Input:
    * --fastq1: the first input FASTQ file
    * --name: the sample name, output files will be named after this name
    * --reference: the reference genome FASTA file, *.fai, *.dict and bwa indexes are required.
    * --gff: the GFFv3 gene annotations file
    * --output: the folder where to publish output

Optional input:
    * --fastq2: the second input FASTQ file
    * --min_base_quality: minimum base call quality to take a base into account (default: 20)
    * --min_mapping_quality: minimum mapping quality to take a read into account (default: 20)
    * --low_frequency_variant_threshold: VAF threshold to mark a variant as low frequency (default: 0.2)
    * --subclonal_variant_threshold: VAF superior threshold to mark a variant as subclonal (default: 0.8)
    * --strand_bias_threshold: threshold for the strand bias test Phred score (default: 20)
    * --memory: the ammount of memory used by each job (default: 3g)
    * --cpus: the number of CPUs used by each job (default: 1)

Output:
    * Output a normalized, phased and annotated VCF file for each of BCFtools, GATK and LoFreq
    * Output a TSV file output from iVar
```


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
