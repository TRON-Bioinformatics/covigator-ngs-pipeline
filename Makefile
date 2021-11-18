
all : clean test

clean:
	rm -rf output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash tests/test_00_help.sh
	#bash tests/test_00_initialize.sh
	bash tests/test_01.sh
	bash tests/test_02.sh
	bash tests/test_03.sh
	bash tests/test_04.sh
	bash tests/test_05.sh
	bash tests/test_06.sh
	bash tests/test_07.sh
	bash tests/test_08.sh
	bash tests/test_09.sh
