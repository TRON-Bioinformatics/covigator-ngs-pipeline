#!/bin/bash

echo "Running CoVigator pipeline test 0 - initialize"
nextflow main.nf -profile test,conda --initialize