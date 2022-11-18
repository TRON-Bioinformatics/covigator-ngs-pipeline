#!/bin/bash

##################################################################################
# FASTA variant caller
##################################################################################
echo "Running CoVigator pipeline Python unit tests"
virtualenv venv --python python3
source venv/bin/activate
pip install -r bin/requirements.txt
python3 -m unittest bin/tests/test_assembly_variant_caller.py
python3 -m unittest bin/tests/test_reference_genome.py
python3 -m unittest bin/tests/test_phasing.py
deactivate
