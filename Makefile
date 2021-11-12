
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash tests/test_00.sh
	bash tests/test_01.sh
	bash tests/test_02.sh
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
	echo "hCoV-19_NTXX\t"`pwd`"/test_data/hCoV-19_NTXX.fasta\n" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --input_fastas_list test_data/test_input.txt \
	--output output/test5
	echo "ERR4145453\t"`pwd`"/test_data/ERR4145453_1.fastq.gz\t"`pwd`"/test_data/ERR4145453_2.fastq.gz\n" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --input_fastqs_list test_data/test_input.txt \
	--library paired --output output/test6
	echo "ERR4145453\t"`pwd`"/test_data/ERR4145453_1.fastq.gz\n" > test_data/test_input.txt
	nextflow main.nf -profile test,conda --input_fastqs_list test_data/test_input.txt \
	--library single --output output/test7

quick_test:
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test1 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz --skip_sarscov2_annotations
	test -s output/test1/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.fastp_stats.json || { echo "Missing test 1 FASTP output file!"; exit 1; }
	test -s output/test1/ERR4145453.fastp_stats.html || { echo "Missing test 1 FASTP output file!"; exit 1; }
	test -s output/test1/ERR4145453.coverage.tsv || { echo "Missing test 1 coverage output file!"; exit 1; }
	test -s output/test1/ERR4145453.depth.tsv || { echo "Missing test 1 depth output file!"; exit 1; }

check:
	test -s output/test1/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.fastp_stats.json || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.fastp_stats.html || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.coverage.tsv || { echo "Missing test 1 coverage output file!"; exit 1; }
	test -s output/test1/ERR4145453.depth.tsv || { echo "Missing test 1 depth output file!"; exit 1; }
	test -s output/test2/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.ivar.tsv || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.fastp_stats.json || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.fastp_stats.html || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.ivar.tsv || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.fastp_stats.json || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453.fastp_stats.html || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test4/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz || { echo "Missing test 4 VCF output file!"; exit 1; }
	test -s output/test5/hCoV-19_NTXX.assembly.normalized.annotated.vcf.gz || { echo "Missing test 5 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.ivar.tsv || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.fastp_stats.json || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test6/ERR4145453.fastp_stats.html || { echo "Missing test 6 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.bcftools.normalized.annotated.vcf.gz || { echo "Missing test 7 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.gatk.normalized.annotated.vcf.gz || { echo "Missing test 7 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.lofreq.normalized.annotated.vcf.gz || { echo "Missing test 7 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.ivar.tsv || { echo "Missing test 7 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.fastp_stats.json || { echo "Missing test 7 VCF output file!"; exit 1; }
	test -s output/test7/ERR4145453.fastp_stats.html || { echo "Missing test 7 VCF output file!"; exit 1; }

test_assembly_variant_caller:
	python3 -m unittest bin/test_assembly_variant_caller.py