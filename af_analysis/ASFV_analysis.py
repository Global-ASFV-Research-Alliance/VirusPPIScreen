import analysis.AFppi as ppi
import analysis.plot as jplot
import pandas as pd
import os, sys, glob, pdb

# protein length databases
ASFV_protlen_df = pd.read_csv("/home/jacob.fenster/scripts/alphafold/data/ASFV-Georgia_2007-LR743116_trans-rename_peptide-len-db.csv", index_col=0)
FMDV_protlen_df = pd.read_csv("/home/jacob.fenster/scripts/alphafold/data/FMDV-AY593768_protein-len-db.csv", index_col=0)
# input directories
FMDV_AFmult_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231220_FMDV_AF-multimer_dimers/"
ASFV_gold_AFmult_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231204_AFtests/ASFV_ASFV-Gold/" 
ASFV_gold_paired_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231220_AF-2.0.1-Ian_paired_unpaired_benchmarks/"
ASFV_gold_paired_unpaired_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231220_AF-2.0.1-Ian_paired_benchmarks/"
ASFV_AFmult_homodimer_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231207_AF_multimer_ASFV_homodimer/"
# output dir
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/"

read_raw_data = False # this is put here to read the raw datasets. Set to false if complete to save time
if read_raw_data:
# extract datasets 
    ASFV_g_mult = ppi.AFmult_data_extract(ASFV_gold_AFmult_dir)
    ASFV_g_ip = ppi.dir_max_contact(ASFV_gold_paired_dir, ASFV_protlen_df)
    ASFV_g_ipu = ppi.dir_max_contact(ASFV_gold_paired_unpaired_dir, ASFV_protlen_df)
    ASFV_homo = ppi.AFmult_data_extract(ASFV_AFmult_homodimer_dir)
    FMDV = ppi.AFmult_data_extract(FMDV_AFmult_dir)
    # merge to get only experimentally validated set of ASFV PPIs
# output raw data
    os.makedirs(output_dir, exist_ok=True), os.makedirs(f"{output_dir}ASFV/", exist_ok=True), os.makedirs(f"{output_dir}FMDV/", exist_ok=True)
    ASFV_g_mult.to_csv(f"{output_dir}/ASFV/ASFV_gold_AFmult_data.csv")
    ASFV_g_ip.to_csv(f"{output_dir}/ASFV/ASFV_gold_Ian2021_paired_data.csv")
    ASFV_g_ipu.to_csv(f"{output_dir}/ASFV/ASFV_gold_Ian2021_paired-unpaired_data.csv")
    ASFV_homo.to_csv(f"{output_dir}/ASFV/ASFV_homodimer_AFmult_data.csv")
    FMDV.to_csv(f"{output_dir}/FMDV/FMDV_AFmult_dimer_data.csv")
else:
# read from saved csv files to save time
    ASFV_g_mult = pd.read_csv(f"{output_dir}/ASFV/ASFV_gold_AFmult_data.csv", index_col=0)
    ASFV_g_ip = pd.read_csv(f"{output_dir}/ASFV/ASFV_gold_Ian2021_paired_data.csv", index_col=0)
    ASFV_g_ipu = pd.read_csv(f"{output_dir}/ASFV/ASFV_gold_Ian2021_paired-unpaired_data.csv", index_col=0)
    ASFV_homo = pd.read_csv(f"{output_dir}/ASFV/ASFV_homodimer_AFmult_data.csv", index_col=0)
    FMDV = pd.read_csv(f"{output_dir}/FMDV/FMDV_AFmult_dimer_data.csv", index_col=0)
    
breakpoint()
ASFV_expg_ip = ppi.custom_merge_any(ASFV_g_mult, ASFV_g_ip, 'PPI', [], ['max_contact_prob'])
ASFV_expg_ipu = ppi.custom_merge_any(ASFV_g_mult, ASFV_g_ipu, 'PPI', [], ['max_contact_prob'])
ASFV_expg_ip.to_csv(f"{output_dir}/ASFV/ASFV_expGold_Ian2021_paired_data.csv")
ASFV_expg_ipu.to_csv(f"{output_dir}/ASFV/ASFV_expGold_Ian2021_paired-unpaired_data.csv")
