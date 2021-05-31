
all : clean test check

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	nextflow main.nf --help
	nextflow main.nf -profile test,conda --name ERR4145453 \
	--fastq1 `pwd`/test_data/ERR4145453_1.fastq.gz \
	--fastq2 `pwd`/test_data/ERR4145453_2.fastq.gz \
	--reference `pwd`/test_data/MN908947.3.fa


check:
	test -s output/ERR4145453/ERR4145453.bam || { echo "Missing BAM output file!"; exit 1; }