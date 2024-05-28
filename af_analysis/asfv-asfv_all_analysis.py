# this script is to analyze the ASFV-ASFV all vs all homo heterodimer set that were processed
# via the Ian et al 2021 AlphaFold script
import os, sys, glob, pdb
import pandas as pd
import analysis.AFppi as ppi
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar

# ASFV Georgia 2007 protein length databases
asfv_protlen_df = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/ASFV/ASFV-Georgia_2007-LR743116_trans-rename_peptide-len-db.csv", index_col=0)
# ASFV-ASFV all vs all AFIan2021 output 
asfv_heterodimer = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/heterodimer/"
asfv_homodimer = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/homodimer/"
# analysis output
project_dir = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
# process AFIan2021 .npz files 
if True:
    asfv_heterodimer_df = ppi.dir_max_contact(asfv_heterodimer, asfv_protlen_df)
    asfv_homodimer_df = ppi.dir_max_contact(asfv_homodimer, asfv_protlen_df)
    asfv_heterodimer_df.to_csv(f"{project_dir}asfv-asfv_heterodimer_maxcontact.csv")
    asfv_homodimer_df.to_csv(f"{project_dir}asfv-asfv_homodimer_maxcontact.csv")



