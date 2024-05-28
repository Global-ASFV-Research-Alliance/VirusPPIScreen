# This script is to analyze the all vs all ASFV, Vaccinia, FMDV data from the mega pMSA pipeline of 
# Ian2021 AlphaFold script. I also will be comparing the AFmultimer results 
# created 20240229

import os, sys, glob, itertools, time, pdb
import pandas as pd
import numpy as np
import analysis.AFppi as ppi
import analysis.pmsa_stats as pstats
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar
import analysis.RFppi as rppi

def protein_len_df(proteome, output_file, uniprot=True):
# Generate a protein length df for max contact analysis 
# fasta formatted proteome file, full path output file, and if it is uniprot format for accession parsing
    descs, seqs = bpar.read_fasta(proteome)
    df = pd.DataFrame(columns=['length'])
    # parse uniprot accession number for naming
    if uniprot:
        descs_new = []
        for desc, seq in zip(descs, seqs):
            name = desc.split('|')[1]
            if seq[-1] == '*':
                length = len(seq) - 1 
            else:
                length = len(seq)
            df.loc[name] = {'length': length}
        df.to_csv(output_file)

def consolodate_output(dirs, output_dir):
    # This is to consolodate Ian2021 output data 
    # dirs is a list of directories that contain output data to consolodate
    for dir in dirs:
        os.system(f"mv {dir}* {output_dir}")

def _load_master(df, master, df_column, master_column):
# loads the df_column value onto the master_column value by finding a PPI_A or PPI_B match
# df.index must be the PPI names in the PPI_A or PPI_B columns of master
    for index in df.index:
        found = False
        for m_index in master.index:
            if index == master.loc[m_index, 'PPI_A'] or index == master.loc[m_index, 'PPI_B']:
                found = True 
                master.loc[m_index, master_column] = df.loc[index, df_column]
        if not found:
            print(f"{index} not found. skipping")
    return master

################## Main ######################
ian2021 = 'Humphreys_maxcontact' # name for Ian Humphrey's 2021 max contact calcuation
afmult = 'iptm+ptm' # name for AFmultimer iptm+ptm calculation
master_tables = False # switch for generating master tables from proteomes etc. 
max_contact = False # switch for calculating max_contact dataframes
afmult_iptm = False # switch for calculating AlphaFold Multimer iptm+ptm scores 
load_master = False # switch for loading parsed data into master tables
neff_analysis = False # switch to do Neffective sequence analysis 
gold_standard = False # add cross reference experimental y2h data etc.  
gold_std_data = False # add gold standard separately run datasets to master dfs
redo_pmsa_stats = False # new neff, avg_pid, avg_cov 
parse_rf = False
load_master_rf = False 
neff_asfv_gold = False 
vacc_homodimer_MINT = True

vacc_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
fmdv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_all_homo-heterodimer/"
asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"

### GENERATE MASTER TABLES. Logic will be that proteinA is always alphabetically before protienB 
columns = ['PPI_A', 'PPI_B', 'proteinA', 'proteinB', ian2021, afmult, 
           'gold_std', 'seqA', 'seqB', 'protA_description', 'protB_description']
if master_tables:
# ASFV master table
    asfv_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/Georgia-2007-LR743116_trans-rename.fa"
    asfv_descriptions_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV-Georgia2007_ProteinDescriptions.csv"

    asfv_proteins, asfv_seqs = bpar.read_fasta(asfv_proteome)
    asfv_prot_desc = pd.read_csv(asfv_descriptions_file, index_col=0)
    for index in asfv_prot_desc.index:
        new_index = index.replace(' ', '_')
        asfv_prot_desc.rename(index={index: new_index}, inplace=True)
    asfv_dict = {} # keys are genes and entries are peptide sequences 
    for desc, seq in zip(asfv_proteins, asfv_seqs):
        asfv_dict[desc] = seq
    # heterodimer 
    asfv_pairs = list(itertools.combinations(asfv_proteins, 2))
    asfv_master = pd.DataFrame(columns=columns)
    for pair, index in zip(asfv_pairs, range(len(asfv_pairs))):
        protA, protB = min(pair), max(pair) 
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': asfv_dict[protA][:-1], 'seqB': asfv_dict[protB][:-1],
                       'protA_description': asfv_prot_desc.loc[protA, 'Description'], 
                       'protB_description': asfv_prot_desc.loc[protB, 'Description']}
        asfv_master.loc[index] = partial_row 
    # homodimer 
    asfv_homodimer_pairs = []
    for protein in asfv_proteins:
        asfv_homodimer_pairs.append((protein, protein))
    for pair, index in zip(asfv_homodimer_pairs, range(len(asfv_master.index), len(asfv_master.index)+len(asfv_homodimer_pairs))): 
        protA, protB = pair[0], pair[1]
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': asfv_dict[protA][:-1], 'seqB': asfv_dict[protB][:-1],
                       'protA_description': asfv_prot_desc.loc[protA, 'Description'], 
                       'protB_description': asfv_prot_desc.loc[protB, 'Description']}
        asfv_master.loc[index] = partial_row 
    # fix MGF_110-10-L_-_MGF110-14L_fusion
    for index in asfv_master.index:
        protA, protB = asfv_master.loc[index, 'proteinA'], asfv_master.loc[index, 'proteinB']
        if protA == 'MGF_110-10-L_-_MGF110-14L_fusion':
            asfv_master.loc[index, 'proteinA'] = 'MGF_110-10-L-MGF110-14L_fusion'
            if protA == protB: 
                asfv_master.loc[index, 'PPI_A'] = 'MGF_110-10-L-MGF110-14L_fusion__MGF_110-10-L-MGF110-14L_fusion'
                asfv_master.loc[index, 'PPI_B'] = 'MGF_110-10-L-MGF110-14L_fusion__MGF_110-10-L-MGF110-14L_fusion'
            else:
                other_protein = asfv_master.loc[index, 'PPI_A'].split('__')[1]
                asfv_master.loc[index, 'PPI_A'] = f'MGF_110-10-L-MGF110-14L_fusion__{other_protein}'
                asfv_master.loc[index, 'PPI_B'] = f'{other_protein}__MGF_110-10-L-MGF110-14L_fusion'
        if protB == 'MGF_110-10-L_-_MGF110-14L_fusion':
            asfv_master.loc[index, 'proteinB'] = 'MGF_110-10-L-MGF110-14L_fusion'
            if protA != protB:
                other_protein = asfv_master.loc[index, 'PPI_A'].split('__')[0]
                asfv_master.loc[index, 'PPI_B'] = f'MGF_110-10-L-MGF110-14L_fusion__{other_protein}'
                asfv_master.loc[index, 'PPI_A'] = f'{other_protein}__MGF_110-10-L-MGF110-14L_fusion'


    asfv_master.to_csv(f"{asfv_analysis_out}asfv_partial_master.csv")
# Vaccinia master table
    vacc_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
    descsi, seqs = bpar.read_fasta(vacc_proteome)
    vacc_dict = {} # Uniprot accession are keys, and secondary keys are seq and Description
    vacc_proteins = [] # just uniprot accessions
    for desc, seq in zip(descsi, seqs):
        accession = desc.split('|')[1]
        vacc_proteins.append(accession)
        vacc_dict[accession] = {'seq': seq, 
                                'Description': " ".join(desc.split('|')[2].split(' OS=')[0].split()[1:])}
    # heterodimers
    vacc_master = pd.DataFrame(columns=columns)
    vacc_pairs = list(itertools.combinations(vacc_proteins, 2))
    for pair, index in zip(vacc_pairs, range(len(vacc_pairs))):
        protA, protB = min(pair), max(pair) 
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': vacc_dict[protA]['seq'], 'seqB': vacc_dict[protB]['seq'], 
                       'protA_description': vacc_dict[protA]['Description'],
                       'protB_description': vacc_dict[protB]['Description']}
        vacc_master.loc[index] = partial_row
    # homodimers
    vacc_homodimer_pairs = []
    for protein in vacc_proteins:
        vacc_homodimer_pairs.append((protein, protein))
    for pair, index in zip(vacc_homodimer_pairs, range(len(vacc_master.index), len(vacc_master.index)+len(vacc_homodimer_pairs))): 
        protA, protB = pair[0], pair[1]
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': vacc_dict[protA]['seq'], 'seqB': vacc_dict[protB]['seq'], 
                       'protA_description': vacc_dict[protA]['Description'],
                       'protB_description': vacc_dict[protB]['Description']}
        vacc_master.loc[index] = partial_row 

    vacc_master.to_csv(f"{vacc_analysis_out}vacc_partial_master.csv")
# FMDV master table 
    fmdv_proteome_file = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/AY593768_proteome.fa"
    descs, seqs = bpar.read_fasta(fmdv_proteome_file)
        # need to find FMDV protein descriptions 
    fmdv_dict = {}
    fmdv_proteins = []
# build fmdv dictionary
    for desc, seq in zip(descs, seqs):
        name = desc.split('68_')[1]
        fmdv_proteins.append(name)
        fmdv_dict[name] = seq
    # heterodimers
    fmdv_master = pd.DataFrame(columns=columns)
    fmdv_pairs = list(itertools.combinations(fmdv_proteins, 2))
    for pair, index in zip(fmdv_pairs, range(len(fmdv_pairs))):
        protA, protB = min(pair), max(pair)
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': fmdv_dict[protA], 'seqB': fmdv_dict[protB]}
        fmdv_master.loc[index] = partial_row
    # homodimers
    fmdv_homodimer_pairs = []
    for protein in fmdv_proteins:
        fmdv_homodimer_pairs.append((protein, protein))
    for pair, index in zip(fmdv_homodimer_pairs, range(len(fmdv_master.index), len(fmdv_master.index)+len(fmdv_proteins))):
        protA, protB = pair[0], pair[1] 
        partial_row = {'PPI_A': f'{protA}__{protB}', 'PPI_B':f'{protB}__{protA}',
                       'proteinA': protA, 'proteinB': protB, 
                       'seqA': fmdv_dict[protA], 'seqB': fmdv_dict[protB]}
        fmdv_master.loc[index] = partial_row
    fmdv_master.to_csv(f"{fmdv_analysis_out}fmdv_partial_master.csv")

### MAX CONTACT from Ian2021 Humphreys and AFmultimer
# Vaccinia analysis 
if max_contact:
    # Generate proteome length df 
    vacc_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
    vacc_lendf_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_protein_len_df.csv"
    protein_len_df(vacc_proteome, vacc_lendf_file, uniprot=True)
    # combine output into a commond folder
    vacc_proj_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/"
    vacc_hetero_dirs = [f"{vacc_proj_out}tranche1/", f"{vacc_proj_out}tranche1_1/", f"{vacc_proj_out}tranche1_2/", 
                            f"{vacc_proj_out}tranche1_3/", f"{vacc_proj_out}tranche1_4/", f"{vacc_proj_out}tranche1_5/", 
                            f"{vacc_proj_out}tranche1_6/", f"{vacc_proj_out}tranche1_6/", f"{vacc_proj_out}tranche2/",
                            f"{vacc_proj_out}tranche3/", f"{vacc_proj_out}tranche4/", f"{vacc_proj_out}tranche5/",
                            f"{vacc_proj_out}tranche6/"]
    vacc_hetero_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/heterodimer/"
    consolodate_output(vacc_hetero_dirs, vacc_hetero_dir)
# process max contact for vaccinia 
    vacc_hetero_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/heterodimer/"
    vacc_homo_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/homodimer/"
    vacc_lendf_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_protein_len_df.csv"

    vacc_len_df = pd.read_csv(vacc_lendf_file, index_col=0)
    vacc_het_df = ppi.dir_max_contact(vacc_hetero_dir, vacc_len_df) 
    vacc_ho_df = ppi.dir_max_contact(vacc_homo_dir, vacc_len_df) 

    vacc_het_df.to_csv(f"{vacc_analysis_out}vaccinia-WR_ian_mega_heterodimer.csv")
    vacc_ho_df.to_csv(f"{vacc_analysis_out}vaccinia-WR_ian_mega_homodimer.csv")
    asfv_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/heterodimer/"
    asfv_ho_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/homodimer/"
    asfv_lendf_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV-Georgia_2007-LR743116_trans-rename_peptide-len-db.csv"

    asfv_len_df = pd.read_csv(asfv_lendf_file, index_col=0)
    asfv_het_df = ppi.dir_max_contact(asfv_het_dir, asfv_len_df)
    asfv_ho_df = ppi.dir_max_contact(asfv_ho_dir, asfv_len_df)
    
    asfv_het_df.to_csv(f"{asfv_analysis_out}asfv-asfv_heterodimer_maxcontact.csv")
    asfv_ho_df.to_csv(f"{asfv_analysis_out}asfv-asfv_homodimer_maxcontact.csv")
# process FMDV data
    fmdv_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/FMDV/Ian2021/heterodimer/"
    fmdv_ho_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/FMDV/Ian2021/homodimer/"
    fmdv_lendf_file = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV-AY593768_protein-len-db.csv"

    fmdv_len_df = pd.read_csv(fmdv_lendf_file, index_col=0)
    fmdv_het_df = ppi.dir_max_contact(fmdv_het_dir, fmdv_len_df)
    fmdv_ho_df = ppi.dir_max_contact(fmdv_ho_dir, fmdv_len_df)

    fmdv_het_df.to_csv(f"{fmdv_analysis_out}fmdv-fmdv_heterodimer_maxcontact.csv")
    fmdv_ho_df.to_csv(f"{fmdv_analysis_out}fmdv-fmdv_homodimer_maxcontact.csv")

### AlphaFold multimer parse high-scoring predictions
if afmult_iptm:
    # parse asfv
    asfv_afmult_hs_het_out  = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/high_score_heterodimer/" 
    asfv_afmult_hs_het = ppi.AFmult_data_extract(asfv_afmult_hs_het_out)
    asfv_afmult_hs_het.to_csv(f"{asfv_analysis_out}asfv_afmult_high-score_heterodimer.csv")
    # parse vaccinia
    vacc_afmult_hs_het_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/afmult_Ian-tophits/"
    vacc_afmult_hs_het = ppi.AFmult_data_extract(vacc_afmult_hs_het_out)
    vacc_afmult_hs_het.to_csv(f"{vacc_analysis_out}vaccinia-WR_afmult_high-score_heterodimer.csv")
    # parse fmdv 
    # This has already be done previously.
    # /lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_all_homo-heterodimer/FMDV_AFmult_all_dimer_data.csv
### Load data
asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_partial_master.csv", index_col=0)
vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_partial_master.csv", index_col=0)
fmdv_master = pd.read_csv(f"{fmdv_analysis_out}fmdv_partial_master.csv", index_col=0)

vacc_het_df = pd.read_csv(f"{vacc_analysis_out}vaccinia-WR_ian_mega_heterodimer.csv", index_col=0) 
vacc_ho_df = pd.read_csv(f"{vacc_analysis_out}vaccinia-WR_ian_mega_homodimer.csv", index_col=0)
asfv_het_df = pd.read_csv(f"{asfv_analysis_out}asfv-asfv_heterodimer_maxcontact.csv", index_col=0)
asfv_ho_df = pd.read_csv(f"{asfv_analysis_out}asfv-asfv_homodimer_maxcontact.csv", index_col=0)
fmdv_het_df = pd.read_csv(f"{fmdv_analysis_out}fmdv-fmdv_heterodimer_maxcontact.csv", index_col=0)
fmdv_ho_df = pd.read_csv(f"{fmdv_analysis_out}fmdv-fmdv_homodimer_maxcontact.csv", index_col=0)

asfv_afmult_hs_het = pd.read_csv(f"{asfv_analysis_out}asfv_afmult_high-score_heterodimer.csv", index_col=0)
asfv_afmult_homodimer = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/ASFV/early_analysis/ASFV_homodimer_AFmult_data.csv", index_col=0)
vacc_afmult_hs_het = pd.read_csv(f"{vacc_analysis_out}vaccinia-WR_afmult_high-score_heterodimer.csv", index_col=0)
fmdv_afmult_all = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_all_homo-heterodimer/FMDV_AFmult_all_dimer_data.csv", index_col=0)

### LOAD DATA into master tables 
if load_master:
# asfv
    timei = time.time()
    print("started loading ian2021 homodimer asfv set")
    # load Ian2021 homodimer data
    asfv_master = _load_master(asfv_ho_df, asfv_master, df_column='max_contact_prob', master_column=ian2021)
    timef = time.time()
    os.system(f"> /home/jacob.fenster/jobs/20240304/{timef-timei:.2f}_elapsed_homodimer.log")
    print(f"finished loading asfv ian2021 ho set. {timef-timei:.2f} seconds  elapsed")
    # load Ian2021 heterodimer data
    print(f"started loading heterodimer asfv set")
    asfv_master = _load_master(asfv_het_df, asfv_master, df_column='max_contact_prob', master_column=ian2021)
    timei = time.time()
    print(f"finished loading asfv heterodimer set. {timei-timef:.2f} seconds elapsed. starting afmult data loading")
    # load alphafold multimer data, heterodimer
    asfv_master = _load_master(asfv_afmult_hs_het, asfv_master, df_column='iptm+ptm', master_column=afmult)
    print('finished afmultimer loading of heterodimer, starting homodimer afmult load')
    # load alphafold multimer data, homodimer
    # clean up names of homodimer set
    for index in asfv_afmult_homodimer.index:
        new_index = f"{index.split('__')[0]}__{index.split('__')[2]}" 
        asfv_afmult_homodimer.rename(index={index:new_index}, inplace=True)
    asfv_master = _load_master(asfv_afmult_homodimer, asfv_master, df_column='iptm+ptm', master_column=afmult)
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial2.csv")
    print('wrote asfv_master_partial2.csv')
# Vaccinia
    # load Ian2021 heterodimer data
    print('starting ian load of vacc master')
    vacc_master = _load_master(vacc_het_df, vacc_master, df_column='max_contact_prob', master_column=ian2021)
    # load Ian2021 homodimer data
    vacc_master = _load_master(vacc_ho_df, vacc_master, df_column='max_contact_prob', master_column=ian2021)
    # load AFmultimer heterodimer data
    vacc_master = _load_master(vacc_afmult_hs_het, vacc_master, df_column='iptm+ptm', master_column=afmult)
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial2.csv")
    print('finished loading vacc data, wrote to output')
# FMDV
    # load Ian2021 heterodimer data
    for index in fmdv_het_df.index:
        new_index = f"{index.split('__')[0].split('768_')[1]}__{index.split('__')[1].split('768_')[1]}"
        fmdv_het_df.rename(index={index:new_index}, inplace=True)
    fmdv_master = _load_master(fmdv_het_df, fmdv_master, df_column='max_contact_prob', master_column=ian2021)
    # load Ian2021 homodimer data 
    for index in fmdv_ho_df.index:
        new_index = f"{index.split('__')[0].split('768_')[1]}__{index.split('__')[1].split('768_')[1]}"
        fmdv_ho_df.rename(index={index:new_index}, inplace=True)
    fmdv_master = _load_master(fmdv_ho_df, fmdv_master, df_column='max_contact_prob', master_column=ian2021)
    # load AFmultimer heterodimer and homodimer data 
    for index in fmdv_afmult_all.index:
        if len(index.split('__')) == 2:
            new_index = f"{index.split('__')[0].split('768_')[1]}__{index.split('__')[1].split('768_')[1]}"
            fmdv_afmult_all.rename(index={index:new_index}, inplace=True)
        elif len(index.split('__')) == 4:
            new_index = f"{index.split('__')[0].split('768_')[1]}__{index.split('__')[2].split('768_')[1]}"
            fmdv_afmult_all.rename(index={index:new_index}, inplace=True)
    fmdv_master = _load_master(fmdv_afmult_all, fmdv_master, df_column='iptm+ptm', master_column=afmult)
    fmdv_master.to_csv(f"{fmdv_analysis_out}fmdv_master_partial2.csv")

# LOAD MASTER DFs 
asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial2.csv", index_col=0)
vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial2.csv", index_col=0)
fmdv_master = pd.read_csv(f"{fmdv_analysis_out}fmdv_master_partial2.csv", index_col=0)

asfv_het_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/pmsas/paired_unpaired/"
vacc_het_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFIan2021/vaccinia-vaccinia-mega/pmsas/paired_unpaired/"
fmdv_het_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/AFIan2021/mega_pMSA/pmsas/paired_unpaired/"

asfv_ho_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/pmsas/homodimer/"
vacc_ho_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFIan2021/vaccinia-vaccinia-mega/pmsas/homodimer/"
fmdv_ho_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/AFIan2021/mega_pMSA/pmsas/homodimer/"

if neff_analysis:
    tmp = "/lustrefs/fadru/projects/asfv-ppi/data/tmp/"
    print("starting neff calc on asfv data")
    asfv_master = pstats.load_neff(asfv_het_pmsa_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)
    asfv_master = pstats.load_neff(asfv_ho_pmsa_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial3.csv")
    print("finished neff calc on asfv data")
    print("starting neff calc on vacc data")
    vacc_master = pstats.load_neff(vacc_het_pmsa_dir, vacc_master, name_switch=None, id_tr=80, tmp=tmp)
    vacc_master = pstats.load_neff(vacc_ho_pmsa_dir, vacc_master, name_switch=None, id_tr=80, tmp=tmp)
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial3.csv")
    print("finished neff calc on vacc data")
    print("starting neff calcl on fmdv data")
    fmdv_master = pstats.load_neff(fmdv_het_pmsa_dir, fmdv_master, name_switch='fmdv', id_tr=80, tmp=tmp)
    fmdv_master = pstats.load_neff(fmdv_ho_pmsa_dir, fmdv_master, name_switch='fmdv', id_tr=80, tmp=tmp)
    fmdv_master.to_csv(f"{fmdv_analysis_out}fmdv_master_partial3.csv")
    print("finished neff calcl on fmdv data")

# Load master df
asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial3.csv", index_col=0)
vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial3.csv", index_col=0)
fmdv_master = pd.read_csv(f"{fmdv_analysis_out}fmdv_master_partial3.csv", index_col=0)

### CROSS REFERENCE WITH EXPERIMENTAL PPI DATA

# make vaccinia McCraith -> uniprot accession dictionary universal
vacc_uniprot_genenames_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_virus-WR_UniProt_genenames.csv"
vacc_mccraith_table_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000_vaccinia_y2h.csv"
vacc_genenames = pd.read_csv(vacc_uniprot_genenames_file, index_col=0)
vacc_mccraith_table = pd.read_csv(vacc_mccraith_table_file)
uniprot_mccraith_cr = {} #cross reference dictionary. Keys are UniProt accession and values are McCraith genenames
mccraith_pairs = set()
for index in vacc_mccraith_table.index:
    protA, protB = vacc_mccraith_table.loc[index, 'protA'], vacc_mccraith_table.loc[index, 'protB']
    mccraith_pairs.add((protA, protB)), mccraith_pairs.add((protB, protA))
    foundA, foundB = False, False
    for accession in vacc_genenames.index:
        if type(vacc_genenames.loc[accession, 'gene_names']) == str:
            for genename in vacc_genenames.loc[accession, 'gene_names'].split(', '):
                if protA == genename:
                    foundA = True
                    uniprot_mccraith_cr[genename] = accession
                if protB == genename:
                    foundB = True
                    uniprot_mccraith_cr[genename] = accession
    if not foundA:
        print(f"Did not find {protA} in uniprot crossref")
    if not foundB:
        print(f"Did not find {protB} in uniprot crossref")

if gold_standard:
    # create McCraith UniProt accession pairs
    mccraith_uniprot_pairs = set()
    for pair in mccraith_pairs:
        try:
           mccraith_uniprot_pairs.add((uniprot_mccraith_cr[pair[0]], uniprot_mccraith_cr[pair[1]]))
        except KeyError:
            print(f"again, {pair} not in uniprot crosref")

    # load vacc master table with McCraith data
    for index in vacc_master.index:
        ppi = vacc_master.loc[index, 'PPI_A']
        pair = (ppi.split('__')[0], ppi.split('__')[1])
        if pair in mccraith_uniprot_pairs:
            vacc_master.loc[index, 'gold_std'] =  1
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial4.csv")

    # asfv exp gold std data
    asfv_gold_std_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV_PPI_gold_stds.csv"
    asfv_gold_stds = pd.read_csv(asfv_gold_std_file)
    asfv_gold_pairs = set()
    for index in asfv_gold_stds.index:
        protA, protB = asfv_gold_stds.loc[index, 'proteinA'], asfv_gold_stds.loc[index, 'proteinB']
        asfv_gold_pairs.add((protA, protB)), asfv_gold_pairs.add((protB, protA))
    for index in asfv_master.index:
        protA, protB = asfv_master.loc[index, 'proteinA'], asfv_master.loc[index, 'proteinB']
        if (protA, protB) in asfv_gold_pairs:
            asfv_master.loc[index, 'gold_std'] = 1
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial4.csv")

    # fmdv y2h work
# NEED TO RUN ADDITIONAL FMDV CLEVAGE VARIANTS 
    fmdv_y2h_file = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_y2h_df.csv"
    fmdv_y2h = pd.read_csv(fmdv_y2h_file, index_col=0)
    for ppi in fmdv_y2h.index:
        found = False
        for index in fmdv_master.index:
            if ppi == fmdv_master.loc[index, 'PPI_A'] or ppi == fmdv_master.loc[index, 'PPI_B']:
                found = True
                fmdv_master.loc[index, 'gold_std'] = 1 
        if not found:
            print(f"{ppi} not found in fmdv_master")
    fmdv_master.to_csv(f"{fmdv_analysis_out}fmdv_master_partial4.csv")

asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial4.csv", index_col=0)
vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial4.csv", index_col=0)
fmdv_master = pd.read_csv(f"{fmdv_analysis_out}fmdv_master_partial4.csv", index_col=0)
# GOLD STD DATASETS
# these were run separately and this will add them to the master dfs
if gold_std_data:
    # asfv 
    asfv_gold_ian_mega_pmsa = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/early_analysis/ASFV_gold_Ian_mega_paired-unpaired.csv"
    asfv_gold_afmult_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/early_analysis/ASFV_gold_AFmult_data.csv"

    asfv_gold_ian = pd.read_csv(asfv_gold_ian_mega_pmsa, index_col=0)
    asfv_gold_afmult = pd.read_csv(asfv_gold_afmult_file, index_col=0)
    asfv_master = _load_master(asfv_gold_ian, asfv_master, df_column='max_contact_prob', master_column=ian2021)
    asfv_master = _load_master(asfv_gold_afmult, asfv_master, df_column='iptm+ptm', master_column=afmult)
    # add Neff data
    asfv_ian_gold_pmsas_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/gold/pmsas/paired_unpaired/"
    tmp = "/lustrefs/fadru/projects/asfv-ppi/data/tmp/"
    print("starting neff calc on asfv data")
    asfv_master = pstats.load_neff(asfv_ian_gold_pmsas_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)

    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial5.csv")
    # vaccinia 
    mccraith_ian_mega_het = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_Ian_heterodimer_paired-unpaired.csv"
    mccraith_ian_mega_ho = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_Ian_homodimer.csv"
    mccraith_afmult_all = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_AFmult.csv"

    mccraith_ian_het = pd.read_csv(mccraith_ian_mega_het, index_col=0)
    mccraith_ian_ho = pd.read_csv(mccraith_ian_mega_ho, index_col=0)
    mccraith_afmult = pd.read_csv(mccraith_afmult_all, index_col=0)

    # the mccraith ian data has already been loaded as it was run twice. keeping files here for unusual names
    # custom load of AFmult where it keeps existing values. This is because the previous run was done with 2x models
    for index in mccraith_afmult.index:
        found = False
        if len(index.split('__')) == 2:
             protA, protB = index.split('__')[0], index.split('__')[1]
        elif len(index.split('__')) == 4:
             protA, protB = index.split('__')[0], index.split('__')[2]
        try:
            ppi = f"{uniprot_mccraith_cr[protA]}__{uniprot_mccraith_cr[protB]}"
            for m_index in vacc_master.index:
                    if ppi == vacc_master.loc[m_index, 'PPI_A'] or ppi == vacc_master.loc[m_index, 'PPI_B']:
                        found = True 
                        if np.isnan(vacc_master.loc[m_index, 'iptm+ptm']):
                            vacc_master.loc[m_index, 'iptm+ptm'] =  mccraith_afmult.loc[index, 'iptm+ptm']
                        else: 
                            print(f"{index} already loaded. Continuing")
        except KeyError:
            print(f"{protA} or {protB} not in uniprot accessions. Skipping")
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial5.csv")

###

asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial5.csv", index_col=0)
vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial5.csv", index_col=0)


tmp = "/lustrefs/fadru/projects/asfv-ppi/data/tmp/"

if redo_pmsa_stats:
    # afmult addition
    # parse asfv
    if False:
        asfv_afmult_hs_het_out  = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/high_score_heterodimer/" 
        asfv_afmult_hs_het = ppi.AFmult_data_extract(asfv_afmult_hs_het_out)
        asfv_afmult_hs_het.to_csv(f"{asfv_analysis_out}asfv_afmult_high-score_heterodimer.csv")
        asfv_hs_het = pd.read_csv(f"{asfv_analysis_out}asfv_afmult_high-score_heterodimer.csv", index_col=0)
        asfv_master = _load_master(asfv_afmult_hs_het, asfv_master, df_column='iptm+ptm', master_column=afmult)
        asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial6.csv")
        asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial6.csv", index_col=0)
### SET INDICIES TO PPI_A
        asfv_master = asfv_master.set_index('PPI_A')
        # asfv pmsa stats
        print('started pMSA stats of asfv het')
        timei = time.time()
        asfv_master = pstats.pmsa_dir_summary_stats(asfv_het_pmsa_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)
        timef = time.time()
        os.system(f'> {asfv_analysis_out}fin_asfv_het_{timef-timei:.1f}.log')
        print('finished pMSA stats of asfv het. starting ASFV homodimer')
        asfv_master = pstats.pmsa_dir_summary_stats(asfv_ho_pmsa_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)
        asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial7.csv")
    vacc_master = vacc_master.set_index('PPI_A')
    # vacc pmsa stats 
    print('started pMSA stats of vacc het')
    vacc_master = pstats.pmsa_dir_summary_stats(vacc_het_pmsa_dir, vacc_master, name_switch=None, id_tr=80, tmp=tmp)
    print('finished pMSA stats of vacc het. staring vacc homodime')
    vacc_master = pstats.pmsa_dir_summary_stats(vacc_ho_pmsa_dir, vacc_master, name_switch=None, id_tr=80, tmp=tmp)

    print('finished pMSA stats of vacc homodimer')
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial6.csv")

### ANALYSIS SWITCHES ###

## local defs ###

def _load_master_index(data_df, master, data_cols, master_cols):
    # Assumes matching indicies (PPI_A)
    # loads data_cols data from data_df onto master_cols master df
    # data_cols and master_cols are lists of column names
    for data_index in data_df.index:
        for data_column, master_column in zip(data_cols, master_cols):
            master.loc[data_index, master_column] = data_df.loc[data_index, data_column]
    return master

def _fix_MGF_fusion(df):
    # changes MGF_110-10-L_-_MGF110-14L_fusion to MGF_110-10-L-MGF110-14L_fusion
    old_mgf_name = 'MGF_110-10-L_-_MGF110-14L_fusion'
    new_mgf_name = 'MGF_110-10-L-MGF110-14L_fusion'
    for index in df.index:
        protA, protB = index.split('__')[0], index.split('__')[1]
        if old_mgf_name == protA and old_mgf_name == protB:
            ppi = f"{new_mgf_name}__{new_mgf_name}" 
            df.rename(index={index: ppi}, inplace=True)
        elif old_mgf_name == protA:
            ppi = f"{new_mgf_name}__{protB}"
            df.rename(index={index: ppi}, inplace=True)
        elif old_mgf_name == protB:
            ppi = f"{protA}__{new_mgf_name}"
            df.rename(index={index: ppi}, inplace=True)
    return df

if parse_rf:
# ASFV RF analysis 
    asfv_rf_mega_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/rosettafold/ASFV-ASFV_mega-pMSA/"
    asfv_rf_mega_het = rppi.rf2t_npz_analysis(asfv_rf_mega_het_dir, symmetry=True, APC=True, name_switch=None)
    asfv_rf_mega_het.to_csv(f"{asfv_analysis_out}rosettafold_asfv_mega-pmsa_heterodimer.csv")
    asfv_rf_paired_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/rosettafold/ASFV-ASFV_clust99/"
    asfv_rf_paired_het = rppi.rf2t_npz_analysis(asfv_rf_paired_het_dir, symmetry=True, APC=True, name_switch=None)
    asfv_rf_paired_het.to_csv(f"{asfv_analysis_out}rosettafold_asfv_paired99_heterodimer.csv")
# Vacc RF analysis 
    vacc_rf_mega_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/rosettafold/vacc-vacc_mega_pmsa_heterodimers/"
    vacc_rf_mega_het = rppi.rf2t_npz_analysis(vacc_rf_mega_het_dir, symmetry=True, APC=True, name_switch=None)
    vacc_rf_mega_het.to_csv(f"{vacc_analysis_out}rosettafold_vacc_mega-pmsa_heterodimer.csv")
    vacc_rf_paired_het_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/rosettafold/vacc-vacc_mega_pmsa_paired_heterodimers/"
    vacc_rf_paired_het = rppi.rf2t_npz_analysis(vacc_rf_paired_het_dir, symmetry=True, APC=True, name_switch=None)
    vacc_rf_paired_het.to_csv(f"{vacc_analysis_out}rosettafold_vacc_paired99_heterodimer.csv")

# load the RoseTTAFold results into the master dfs
if load_master_rf:
    asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial7.csv", index_col=0)
    vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial6.csv", index_col=0)
    # load RF data
    asfv_rf_mega_het = pd.read_csv(f"{asfv_analysis_out}rosettafold_asfv_mega-pmsa_heterodimer.csv", index_col=0)
    asfv_rf_paired_het = pd.read_csv(f"{asfv_analysis_out}rosettafold_asfv_paired99_heterodimer.csv", index_col=0)
    # fix MGF fusion naming
    asfv_rf_mega_het = _fix_MGF_fusion(asfv_rf_mega_het)
    asfv_rf_paired_het = _fix_MGF_fusion(asfv_rf_paired_het)

    vacc_rf_mega_het = pd.read_csv(f"{vacc_analysis_out}rosettafold_vacc_mega-pmsa_heterodimer.csv", index_col=0)
    vacc_rf_paired_het = pd.read_csv(f"{vacc_analysis_out}rosettafold_vacc_paired99_heterodimer.csv", index_col=0)

    data_cols = ['RF2t++ scores']
    # load asfv
    print(f"started loading asfv df")
    asfv_master = _load_master_index(asfv_rf_mega_het, asfv_master, data_cols, ['RF2++_mega-pMSA'])
    asfv_master = _load_master_index(asfv_rf_paired_het, asfv_master, data_cols, ['RF2++_paired-pMSA'])
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial8.csv")
    print(f"finished loading asfv df")
    print(f'started loading vacc df')
    vacc_master = _load_master_index(vacc_rf_mega_het, vacc_master, data_cols, ['RF2++_mega-pMSA'])
    vacc_master = _load_master_index(vacc_rf_paired_het, vacc_master, data_cols, ['RF2++_paired-pMSA'])
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial7.csv")
    print('finished loading vacc df')
    
if neff_asfv_gold: 
    # partial9 added 'homodimer' column which is 0 for heterodimer and 1 for homodimer 
    asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial9.csv", index_col=0)
    asfv_gold_mega_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/gold/pmsas/paired_unpaired/"
    asfv_master = pstats.pmsa_dir_summary_stats(asfv_gold_mega_pmsa_dir, asfv_master, name_switch=None, id_tr=80, tmp=tmp)
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial9.csv")

# load the vacc homodimer high score and MINT positives 
if vacc_homodimer_MINT:
    vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial9.csv", index_col=0)

    vacc_afmult_hs_ho_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/afmult_ian_top_homodimers/" 
    vacc_afmult_mint_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/MINT_positives/AFmult/"
    vacc_afmult_hs_ho = ppi.AFmult_data_extract(vacc_afmult_hs_ho_dir)
    vacc_afmult_mint = ppi.AFmult_data_extract(vacc_afmult_mint_dir)
    breakpoint()
    vacc_master = _load_master_index(vacc_afmult_hs_ho, vacc_master, ['iptm+ptm'], [afmult])
    vacc_master = _load_master_index(vacc_afmult_mint, vacc_master, ['iptm+ptm'], [afmult])
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial10.csv")
