
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test1 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test2 \
	--fastq1 test_data/ERR4145453_1.fastq.gz
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output output/test3 \
	--fastq1 test_data/ERR4145453_1.fastq.gz \
	--fastq2 test_data/ERR4145453_2.fastq.gz \
	--keep_intermediate

check:
	test -s output/test1/ERR4145453/ERR4145453.bcftools.normalized.phased.annotated.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.gatk.normalized.phased.annotated.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.lofreq.normalized.phased.annotated.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.bcftools.normalized.phased.annotated.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.gatk.normalized.phased.annotated.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.lofreq.normalized.phased.annotated.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.bcftools.normalized.phased.annotated.vcf || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.gatk.normalized.phased.annotated.vcf || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.lofreq.normalized.phased.annotated.vcf || { echo "Missing test 3 VCF output file!"; exit 1; }
	test -s output/test3/ERR4145453/ERR4145453.ivar.tsv || { echo "Missing test 3 VCF output file!"; exit 1; }