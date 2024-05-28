import analysis.AFppi as ppi
import analysis.plot as jplot
import pandas as pd
import os, sys, glob, pdb

# protein length databases
asfv_protlen_df = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV-Georgia_2007-LR743116_trans-rename_peptide-len-db.csv", index_col=0)
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/ASFV/"
asfv_g_mega_pu = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/gold/AF-Ian2021/mega_pMSA_paired_unpaired/"

asfv_g_mega_pu_df = ppi.dir_max_contact(asfv_g_mega_pu, asfv_protlen_df)
asfv_g_mega_pu_df.to_csv(f"{output_dir}ASFV_gold_Ian_mega_paired-unpaired.csv")
