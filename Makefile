
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --initialize
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test1 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test2 \
	--fastq1 test_data/ERR4145453_1.fastq.gz --keep_intermediate
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test3 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz \
	--keep_intermediate
	nextflow main.nf -profile test,conda --name hCoV-19_NTXX \
	--output output/test4 \
	--fasta test_data/hCoV-19_NTXX.fasta
	#python3 -m unittest bin/test_assembly_variant_caller.py

quick_test:
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test1 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz --skip_sarscov2_annotations
	test -s output/test1/ERR4145453/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.fastp_stats.json || { echo "Missing test 1 FASTP output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.fastp_stats.html || { echo "Missing test 1 FASTP output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.coverage.tsv || { echo "Missing test 1 coverage output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.depth.tsv || { echo "Missing test 1 depth output file!"; exit 1; }

check:
	test -s output/test1/ERR4145453/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.fastp_stats.json || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.fastp_stats.html || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.coverage.tsv || { echo "Missing test 1 coverage output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.depth.tsv || { echo "Missing test 1 depth output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.fastp_stats.json || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.fastp_stats.html || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.fastp_stats.json || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.fastp_stats.html || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test4/hCoV-19_NTXX/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz || { echo "Missing test 4 VCF output file!"; exit 1; }