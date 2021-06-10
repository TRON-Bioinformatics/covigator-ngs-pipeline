# Covigator NGS pipeline

The Covigator NGS pipeline process SARS-CoV-2 FASTQ files into analysis ready VCF files. The pipeline is implemented in the Nextflow framework (Di Tommaso, 2017).

The pipeline includes the following steps:
- **Alignment**. `BWA mem` is used for the alignment of single or paired end samples.
- **BAM preprocessing**. BAM files are prepared and duplicate reads are marked using GATK and Picard tools.
- **Variant calling**. `bcftools mpileup` and `bcftools call` is employed to call variants.
- **Variant normalization**. `bcftools norm` and `vt` tools are employed to left align indels, trim variant calls and remove variant duplicates.
- **Variant consequence annotation**. `bcftools csq` is employed to annotate the variant consequences of variants.

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
    * --memory: the ammount of memory used by each job (default: 3g)
    * --cpus: the number of CPUs used by each job (default: 1)

Output:
    * Output a normalized and annotated VCF file.
```


## References

- Di Tommaso, P., Chatzou, M., Floden, E. W., Barja, P. P., Palumbo, E., & Notredame, C. (2017). Nextflow enables reproducible computational workflows. Nature Biotechnology, 35(4), 316–319. https://doi.org/10.1038/nbt.3820
- Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60. [PMID: 19451168]
- Adrian Tan, Gonçalo R. Abecasis and Hyun Min Kang. Unified Representation of Genetic Variants. Bioinformatics (2015) 31(13): 2202-2204](http://bioinformatics.oxfordjournals.org/content/31/13/2202) and uses bcftools [Li, H. (2011). A statistical framework for SNP calling, mutation discovery, association mapping and population genetical parameter estimation from sequencing data. Bioinformatics (Oxford, England), 27(21), 2987–2993. 10.1093/bioinformatics/btr509
- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H. Twelve years of SAMtools and BCFtools. Gigascience. 2021 Feb 16;10(2):giab008. doi: 10.1093/gigascience/giab008. PMID: 33590861; PMCID: PMC7931819.
- Van der Auwera GA, Carneiro M, Hartl C, Poplin R, del Angel G, Levy-Moonshine A, Jordan T, Shakir K, Roazen D, Thibault J, Banks E, Garimella K, Altshuler D, Gabriel S, DePristo M. (2013). From FastQ Data to High-Confidence Variant Calls: The Genome Analysis Toolkit Best Practices Pipeline. Curr Protoc Bioinformatics, 43:11.10.1-11.10.33. DOI: 10.1002/0471250953.bi1110s43.
