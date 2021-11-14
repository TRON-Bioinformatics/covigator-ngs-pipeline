#!/bin/bash

##################################################################################
# FASTA variant caller
##################################################################################
echo "Running CoVigator pipeline test 9"
virtualenv venv --python python3
source venv/bin/activate
pip install biopython==1.79
python3 -m unittest bin/test_assembly_variant_caller.py
deactivate
