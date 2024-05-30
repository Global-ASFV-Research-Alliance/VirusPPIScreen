#!/bin/bash
# This script is to test the pMSA generation pipeline with test files in the examples/example_data/ directory
# Usage:
# cd ViursPPIScreen/ 
# chmod +x ./examples/test_pmsa.sh
# ./examples/test_pmsa.sh
        
command="python3 pmsa_gen.py \
--proteome_file examples/example_data/ASFV_test_genome.fa \
--proteome_type raw \
--strain_diverse_db examples/example_data/example_viral_db.fasta \
--alphafold_msa_dir examples/example_data/AlphaFold_MSAs/ \
--hhsuite_dir hhsuite/ \
--pmsa_output_dir examples/example_pmsa_out/ \
--cpus 2"
echo $command
$command
