# this is to generate the FMDV-FMDV homo-heterodimer set using the mega pMSA algorithm

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
    combined_msa_out = f"{project_dir}combined_msas/"
    temp = f"{project_dir}tmp/"
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
        if False:
            with open(f"{query_fasta_dir}{basename}.fasta", 'w') as f:
                f.write(f">{protein}\n{seq}")
    if False:
        query_fastas = glob.glob(f"{query_fasta_dir}*.fasta")
        samples = len(query_fastas)
        all_NCBI_dbs = repeat(all_NCBI_db, samples)
        jhmmer_outs = repeat(jhmmer_out, samples)
        temps = repeat(temp, samples)
        jhmmer_threads = 8
        jhmmer_threadss = repeat(jhmmer_threads, samples)
        iterations = 3
        iterationss = repeat(iterations, samples)
        jhmmers = repeat('jackhmmer', samples)
        num_workers = 12
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            jhmmer_results = executor.map(pmsa.sub_jackhmmer_search, query_fastas, all_NCBI_dbs, jhmmer_outs, 
                                          temps, jhmmer_threadss, iterationss, jhmmers)
# generate combined MSAs
    if False:
        samples = len(protein_names)
        combined_msa_outs = repeat(combined_msa_out, samples)
        jhmmer_dirs = repeat(jhmmer_out, samples)
        query_fasta_dirs = repeat(query_fasta_dir, samples)
        af_outs = repeat(afmult_out, samples)
        temps = repeat(temp, samples)
        reformats = repeat(reformat, samples)
        num_workers = 12
        #pmsa.combined_msa(protein_names[0], combined_msa_out, jhmmer_out, query_fasta_dir, vaccinia_homodimer_af_out, temp, reformat)
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(pmsa.combined_msa, protein_names, combined_msa_outs, jhmmer_dirs, 
                                       query_fasta_dirs, af_outs, temps, reformats)
# generate pMSAs
    if True:
        # heterodimers 
        pairs = list(itertools.combinations(protein_names, 2)) #all pairwise combinations of input files
        proteinA_msas, proteinB_msas = [], []
        for pair in pairs:
            proteinA_msa = f"{combined_msa_out}{pair[0]}.a3m"
            proteinB_msa = f"{combined_msa_out}{pair[1]}.a3m"
            pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99,
                              unpaired_tr=90, homodimer_tr=95, temp=temp, hhfilter=hhfilter)

        # generate homodimers
        dimers = []
        proteinA_msas, proteinB_msas = [], []
        for protein in protein_names:
            dimers.append((protein, protein))
        for dimer in dimers:
            proteinA_msa = f"{combined_msa_out}{dimer[0]}.a3m"
            proteinB_msa = f"{combined_msa_out}{dimer[0]}.a3m"
            pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99,
                              unpaired_tr=90, homodimer_tr=95, temp=temp, hhfilter=hhfilter)

