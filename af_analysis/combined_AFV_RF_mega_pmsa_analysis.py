# This script is to analyze the RoseTTAFold data associated with the ASFV, vaccinia, and FMDV work
import os, sys, glob, pdb
import numpy as np
import pandas as pd
import analysis.RFppi as rppi

### ANALYSIS SWITCHES ###
parse_rf = False
load_master = True ##MUST RENAME MASTER TABLES

vacc_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
fmdv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_all_homo-heterodimer/"
asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
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


if load_master:
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
    if False:
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
    



    
