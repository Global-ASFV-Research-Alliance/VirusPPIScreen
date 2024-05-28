# This script generates an all vs all ASFV-ASFV heterodimer and homodimer pMSA 
# The algorithm generates a mega combined pMSA from AF output 
import os, glob, pdb, itertools
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat

def _combined_msa(protein_name, combined_msa_out, jhmmer_dir, query_fasta_dir, afmult_out, temp, reformat):
    # generates a combined MSA. Def to be parallelized 
    print(f"starting combined MSA for {protein_name}")
    jhmmer_file = f"{jhmmer_dir}{protein_name}_jhmmer.sto"
    query_fasta = f"{query_fasta_dir}{protein_name}.fasta"
    master_seqs, master_descs = [], []
    viral_descs, viral_seqs = pmsa.cleanup_sto_or_fasta(jhmmer_file, gap_ratio_tr=0.5, genome=True)
    print(f"read jackhmmer sto file")
    if viral_descs[0] == 'query': # fix old jhmmer files that had 'query' as query
        viral_descs[0] = protein_name
    viral_seqs = pmsa.un_align_sto_or_fasta(viral_seqs) 
    master_descs.extend(viral_descs), master_seqs.extend(viral_seqs)
    # read and combine  MSAs from AFmult output
    af_descs, af_seqs = pmsa.combine_AF_MSAs(protein_name, afmult_out, temp=temp, reformat=reformat)
    print("read AF MSAs")
    if af_descs is None:
        print(f"No AFmult msas for {protein}. skipping...")
    else:
        master_descs.extend(af_descs), master_seqs.extend(af_seqs)
        aln_descs, aln_seqs = pmsa.sys_gen_seed_hmmalign(query_fasta, master_descs, master_seqs, temp=f"{temp}aln/", phmmer='phmmer', hmmbuild='hmmbuild', hmmalign='hmmalign', reformat=reformat)
        print("generated alignment")
        # remove duplicate alignments prioritizing those in the all NCBI virus db
        unique_aln_descs, unique_aln_seqs = pmsa.remove_duplicates_prioitize_NCBI(aln_descs, aln_seqs)
        # send a3m MSA to output file  
        bpar.fasta_out_single_line(unique_aln_descs, unique_aln_seqs, f"{combined_msa_out}{protein_name}.a3m")
        print(f"Combined MSA for {protein_name} complete")

# generate homodimers 
if __name__ == '__main__':
# generate output directories
    project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/"
    query_fasta_dir = f"{project_dir}query_fastas/"
    combined_msa_out = f"{project_dir}combined_msas/"
    temp = f"{project_dir}tmp/"
    pmsa_out = f"{project_dir}pmsas/"
    os.makedirs(project_dir, exist_ok=True), os.makedirs(query_fasta_dir, exist_ok=True)
    os.makedirs(combined_msa_out, exist_ok=True), os.makedirs(temp, exist_ok=True)
    os.makedirs(pmsa_out, exist_ok=True)
# software input section
    reformat = "/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl"
    hhfilter = "/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter"

# input section
    jhmmer_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/jhmmer_all_NCBI/"
    asfv_homodimer_out = "/lustrefs/fadru/projects/asfv-ppi/Output/20231207_AF_multimer_ASFV_homodimer/"
    asfv_homodimer_fasta_in = "/lustrefs/fadru/projects/asfv-ppi/Input/20231207_AF_multimer_homodimer/"

# exclude large proteins that we are not interested in 
    exclude_set = {'CP2475L', 'NP1450L', 'G1340L', 'M1249L', 'EP1242L', 'G1211R', 'P1192R',\
            'D1133L', 'F1055L', 'C962R', 'B962L', 'A859L'}
# generate query fastas and protein set
    asfv_proteins = [] 
    asfv_fastas = glob.glob(f"{asfv_homodimer_fasta_in}*.fasta")
    for fasta in asfv_fastas:
        descs, seqs = bpar.read_fasta(fasta)
        protein = descs[0].split('__')[0]
        if protein not in exclude_set:
            asfv_proteins.append(protein)
            if False:
                bpar.fasta_out_single_line([protein], [seqs[0]], f"{query_fasta_dir}{protein}.fasta")
# find list of MSAs
# generate combined MSAs using parallel computing 
    # set number of cores
    num_workers = 2 
    # set constants to iteratables 
    samples = len(asfv_proteins)
    combined_msa_outs, jhmmer_dirs, query_fasta_dirs, temps, reformats = repeat(combined_msa_out, samples),\
            repeat(jhmmer_dir, samples), repeat(query_fasta_dir, samples), repeat(temp, samples), repeat(reformat, samples)
    asfv_homodimer_outs = repeat(asfv_homodimer_out, samples)
    # compute combined MSAs in parallel
    print(f"Starting combined MSA generation for {samples} MSAs")
    if False:
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(_combined_msa, asfv_proteins, combined_msa_outs, jhmmer_dirs, query_fasta_dirs, asfv_homodimer_outs, temps, reformats)
        print("Combined MSA generation complete")
# trying to re-run the failed proteins
# compute pMSAs in parallel 
    # generate hetero dimers
    pairs = list(itertools.combinations(asfv_proteins, 2)) #all pairwise combinations of input files
    proteinA_msas, proteinB_msas = [], []
    for pair in pairs:
        proteinA_msas.append(f"{combined_msa_out}{pair[0]}.a3m")
        proteinB_msas.append(f"{combined_msa_out}{pair[1]}.a3m")
    # make constants iterables
    samples = len(pairs)
    output_dirs, temps, hhfilters = repeat(pmsa_out, samples), repeat(temp, samples), repeat(hhfilter, samples)
    paired_trs = repeat(99, samples)
    unpaired_trs = repeat(90, samples) 
    homodimer_trs = repeat(95, samples)
    if True: #wait to start pMSA generation
        print(f"starting pMSA generation of {samples} heteromers")
        if False:
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                results = executor.map(pmsa.generate_pMSA, proteinA_msas, proteinB_msas, output_dirs, paired_trs, unpaired_trs, homodimer_trs, temps, hhfilters)
        # conduct linear pMSA gen to get error messages 
        if False:
            os.makedirs(f"{project_dir}pmsa_sleuth", exist_ok=True), os.makedirs(f"{project_dir}tmp2/", exist_ok=True)
            for proteinA_msa, proteinB_msa in zip(proteinA_msas, proteinB_msas):
                pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, f"{project_dir}pmsa_sleuth/", 99, 90, 95, f"{project_dir}tmp2/", hhfilter)
                print('')


        print("pMSA generation of heteromers complete")
        # generate homodimers
        dimers = []
        proteinA_msas, proteinB_msas = [], []
        for protein in asfv_proteins:
            dimers.append((protein, protein))
        for dimer in dimers:
            proteinA_msas.append(f"{combined_msa_out}{dimer[0]}.a3m")
            proteinB_msas.append(f"{combined_msa_out}{dimer[0]}.a3m")
        # make constants iterables
        samples = len(dimers)
        output_dirs, temps, hhfilters = repeat(pmsa_out, samples), repeat(temp, samples), repeat(hhfilter, samples)
        paired_trs = repeat(99, samples)
        unpaired_trs = repeat(90, samples) 
        homodimer_trs = repeat(95, samples)
        print(f"starting generation of {samples} homomers")
        if False:
            with ProcessPoolExecutor(max_workers=num_workers) as executor:
                results = executor.map(pmsa.generate_pMSA, proteinA_msas, proteinB_msas, output_dirs, paired_trs, unpaired_trs, homodimer_trs, temps, hhfilters)
        if True:
            for proteinA_msa, proteinB_msa in zip(proteinA_msas, proteinB_msas):
                pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, f"{project_dir}pmsas/", 99, 90, 95, f"{project_dir}tmp2/", hhfilter)
        print(f"ASFV-ASFV_mega_pmsa.py complete")    

