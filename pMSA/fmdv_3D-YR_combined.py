# This is to run just the #3D-YR protein to see if more CPU access increases alignment time
import os, glob, pdb, itertools
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat


if __name__ == '__main__':
    all_NCBI_db = "/lustrefs/fadru/projects/asfv-ppi/data/CompleteGenome_virus_CDS_NCBI_v2_cleanup.fasta"
# output directories  
    project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/AFIan2021/mega_pMSA/" 
    query_fasta_dir = f"{project_dir}query_fasta/"
    jhmmer_out = f"{project_dir}jackhmmer_all_virus/"
    combined_msa_out = f"{project_dir}combined_msas2/"
    temp = f"{project_dir}tmp2/"
    pmsa_out = f"{project_dir}pmsas/"
    os.makedirs(project_dir, exist_ok=True), os.makedirs(query_fasta_dir, exist_ok=True), os.makedirs(jhmmer_out, exist_ok=True)
    os.makedirs(combined_msa_out, exist_ok=True), os.makedirs(temp, exist_ok=True), os.makedirs(pmsa_out, exist_ok=True)

# hmmer and hh-suite software input section
    reformat = "/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl"
    hhfilter = "/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter"

    fmdv_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/AY593768_proteome_mulitimer_input.fa"
    afmult_out = "/lustrefs/fadru/projects/asfv-ppi/Output/FMDV/alphafold-multimer/hetero_homodimers/"

    protein_names = []
    descs, seqs = bpar.read_fasta(fmdv_proteome)
    for desc, seq in zip(descs, seqs):
        protein, basename = desc, desc
        protein_names.append(protein)

    if True:
        protein_name = 'AY593768_3D_YR'
        pmsa.combined_msa(protein_name, combined_msa_out, jhmmer_out, query_fasta_dir, afmult_out, temp, reformat)
