# This script takes outputs of AF and combines them with the NCBI all virus database jackhmmer search
# to generate paired and paired-unpaired pMSAs with maximal diversity 
# Process
import os, glob, pdb
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa

# Generate query fasta sequences for the vaccinia dataset 
asfv_gold_afmult_out_parent = "/lustrefs/fadru/projects/asfv-ppi/Output/20231204_AFtests/ASFV_ASFV-Gold/"
# Input databases, gold standards, jhmmer runs, and proteome data 
all_NCBI_db = "/lustrefs/fadru/projects/asfv-ppi/data/CompleteGenome_virus_CDS_NCBI_v2_cleanup.fasta"
jhmmer_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/jhmmer_all_NCBI/"
gold_csv = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV_PPI_gold_stds.csv"
asfv_proteome_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/Georgia-2007-LR743116_trans-rename.fa"
# output directories  
project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/gold/"
query_fasta_dir = f"{project_dir}query_fastas/"
combined_msa_out = f"{project_dir}combined_msas/"
temp = f"{project_dir}tmp/"
pmsa_out = f"{project_dir}pmsas/"

os.makedirs(project_dir, exist_ok=True), os.makedirs(query_fasta_dir, exist_ok=True)
os.makedirs(combined_msa_out, exist_ok=True), os.makedirs(temp, exist_ok=True)
os.makedirs(pmsa_out, exist_ok=True)
# hmmer and hh-suite software input section
reformat = "/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl"
hhfilter = "/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter"

pairs = [] # holds pairs from .csv input
# generate query fasta files and gold pairs 
gold_df = pd.read_csv(gold_csv)
asfv_descs, asfv_seqs = bpar.read_fasta(asfv_proteome_file) 
asfv_proteome = {}
for desc, seq in zip(asfv_descs, asfv_seqs):
    asfv_proteome[desc] = seq[:-1]
for index in gold_df.index:
    proteinA, proteinB = gold_df.loc[index, 'proteinA'], gold_df.loc[index, 'proteinB']
    pairs.append((proteinA, proteinB))
    with open(f"{query_fasta_dir}{proteinA}.fasta", 'w') as f:
        try:
            f.write(f">{proteinA}\n{asfv_proteome[proteinA]}\n")
        except:
            breakpoint()
    with open(f"{query_fasta_dir}{proteinB}.fasta", 'w') as f:
        try:
            f.write(f">{proteinB}\n{asfv_proteome[proteinB]}\n")
        except:
            breakpoint()

# generate combined msas
for query_fasta in glob.glob(f"{query_fasta_dir}*.fasta"):
    protein_name = os.path.basename(query_fasta).split('.')[0]
    jhmmer_file = f"{jhmmer_dir}{protein_name}_jhmmer.sto"
    master_seqs, master_descs = [], []
    viral_descs, viral_seqs = pmsa.cleanup_sto_or_fasta(jhmmer_file, gap_ratio_tr=0.5, genome=True)
    if viral_descs[0] == 'query': # fix old jhmmer files that had 'query' as query
        viral_descs[0] = protein_name
    viral_seqs = pmsa.un_align_sto_or_fasta(viral_seqs) 
    master_descs.extend(viral_descs), master_seqs.extend(viral_seqs)
    # read and combine  MSAs from AFmult output
    af_descs, af_seqs = pmsa.combine_AF_MSAs(protein_name, asfv_gold_afmult_out_parent, temp=temp, reformat=reformat)
    if af_descs is None:
        print(f"No AFmult msas for {protein}. skipping...")
        breakpoint()
        continue
    else:
        master_descs.extend(af_descs), master_seqs.extend(af_seqs)
        aln_descs, aln_seqs = pmsa.sys_gen_seed_hmmalign(query_fasta, master_descs, master_seqs, temp=f"{temp}aln/", phmmer='phmmer', hmmbuild='hmmbuild', hmmalign='hmmalign', reformat=reformat)
        # remove duplicate alignments prioritizing those in the all NCBI virus db
        unique_aln_descs, unique_aln_seqs = pmsa.remove_duplicates_prioitize_NCBI(aln_descs, aln_seqs)
        # send a3m MSA to output file  
        bpar.fasta_out_single_line(unique_aln_descs, unique_aln_seqs, f"{combined_msa_out}{protein_name}.a3m")

for pair in pairs:
    proteinA_msa = f"{combined_msa_out}{pair[0]}.a3m"
    proteinB_msa = f"{combined_msa_out}{pair[1]}.a3m"
    # generate paired and unpaired pMSAs and ouptut to .a3m files. This returns the raw data for downstream stats
    pmsa.generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99, unpaired_tr=90, homodimer_tr=95, temp=temp,  hhfilter=hhfilter)


