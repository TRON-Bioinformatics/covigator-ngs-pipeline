
all : clean test

clean:
	rm -rf tests/output
	rm -f .nextflow.log*
	rm -rf .nextflow*

test:
	bash tests/scripts/test_00_help.sh
	#bash tests/scripts/test_00_initialize.sh
	bash tests/scripts/test_01.sh
	bash tests/scripts/test_02.sh
	bash tests/scripts/test_03.sh
	bash tests/scripts/test_04.sh
	bash tests/scripts/test_05.sh
	bash tests/scripts/test_06.sh
	bash tests/scripts/test_07.sh
	bash tests/scripts/test_08.sh
	bash tests/scripts/test_09.sh
	bash tests/scripts/test_10.sh
	bash tests/scripts/test_python_unit_tests.sh
