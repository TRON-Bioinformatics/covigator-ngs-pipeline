
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output `pwd`/output/test1 \
	--fastq1 `pwd`/test_data/ERR4145453_1.fastq.gz \
	--fastq2 `pwd`/test_data/ERR4145453_2.fastq.gz
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--output `pwd`/output/test2 \
	--fastq1 `pwd`/test_data/ERR4145453_1.fastq.gz


check:
	test -s output/test1/ERR4145453.bam || { echo "Missing test 1 BAM output file!"; exit 1; }
	test -s output/test1/ERR4145453.preprocessed.bam || { echo "Missing test 1 BAM output file!"; exit 1; }
	test -s output/test1/ERR4145453.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.normalized.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test1/ERR4145453.annotated.vcf || { echo "Missing test 1 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.bam || { echo "Missing test 2 BAM output file!"; exit 1; }
	test -s output/test2/ERR4145453.preprocessed.bam || { echo "Missing test 2 BAM output file!"; exit 1; }
	test -s output/test2/ERR4145453.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.normalized.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }
	test -s output/test2/ERR4145453.annotated.vcf || { echo "Missing test 2 VCF output file!"; exit 1; }