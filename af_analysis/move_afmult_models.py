# This script is to move the pdb models or other files of a given PPI from an AFmultimer output
import os, sys, pdb, glob
import pandas as pd

afmult_ho_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/homodimer/"
afmult_het_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/high_score_heterodimer/"

asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
filtered_files_out = f"{asfv_analysis_out}ASFV_all_AFmult_models/"
asfv_master = pd.read_csv(f"/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/asfv_master_partial10.csv", index_col=0)
afmult_dirs = glob.glob(f"{afmult_het_output_dir}*/")
afmult_dirs = afmult_dirs + glob.glob(f"{afmult_ho_output_dir}*/")
tbs_df = pd.DataFrame(columns=asfv_master.columns)
ppis = set(asfv_master.index.tolist())


if True:
    os.makedirs(filtered_files_out, exist_ok=True)
    for dir in afmult_dirs:
        prot1, prot2 = os.path.dirname(dir).split('/')[-1].split('__')[0], os.path.dirname(dir).split('/')[-1].split('__')[1]
        protA, protB = min([prot1, prot2]), max([prot1, prot2])
        ppi = f'{protA}__{protB}'
        try:
            os.system(f"cp {dir}ranked_0.pdb {filtered_files_out}{ppi}-{asfv_master.loc[ppi,'iptm+ptm']*100:.0f}-ranked_0.pdb")
        except:
            breakpoint()

if False:
    afmult_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/homodimer/"

    asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
    filtered_files_out = f"{asfv_analysis_out}top_hits_homodimer/"
    asfv_master = pd.read_csv(f"/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/asfv_master_partial9.csv", index_col=0)
    afmult_dirs = glob.glob(f"{afmult_output_dir}*/")
    os.makedirs(filtered_files_out, exist_ok=True)

    high_score_het = asfv_master[(asfv_master['iptm+ptm'] > 0.75) & (asfv_master['homodimer']==0)]
    high_score_ho = asfv_master[(asfv_master['iptm+ptm'] > 0.75) & (asfv_master['homodimer']==1)]
    for ppi in high_score_ho.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1]
        afmult_score = high_score_ho.loc[ppi, 'iptm+ptm']
        for dir in afmult_dirs:
            prot1, prot2 = os.path.dirname(dir).split('/')[-1].split('__')[0], os.path.dirname(dir).split('/')[-1].split('__')[1]
            protA_dir, protB_dir = min([prot1, prot2]), max([prot1, prot2])
            if f"{protA}__{protB}" == f"{protA_dir}__{protB_dir}":
                os.system(f"cp {dir}ranked_0.pdb {filtered_files_out}{prot1}__{prot2}-{afmult_score*100:.0f}-ranked_0.pdb")
# vaccinia high score afmult
    afmult_het_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/afmult_Ian-tophits/"
    afmult_ho_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/afmult_ian_top_homodimers/"
    vacc_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"

    filtered_files_out = f"{vacc_analysis_out}Vaccinia_all_AFmult_models/"
    vacc_master = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/vacc_master_partial10.csv", index_col=0)
    afmult_dirs = glob.glob(f"{afmult_het_output_dir}*/")
    afmult_dirs = afmult_dirs + glob.glob(f"{afmult_ho_output_dir}*/")
    os.makedirs(filtered_files_out, exist_ok=True)
    # all pdb files
    for dir in afmult_dirs:
        prot1, prot2 = os.path.dirname(dir).split('/')[-1].split('__')[0], os.path.dirname(dir).split('/')[-1].split('__')[1]
        protA, protB = min([prot1, prot2]), max([prot1, prot2])
        ppi = f'{protA}__{protB}'
        os.system(f"cp {dir}ranked_0.pdb {filtered_files_out}{ppi}-{vacc_master.loc[ppi,'iptm+ptm']*100:.0f}-ranked_0.pdb")

if False:
    high_score = vacc_master[vacc_master['iptm+ptm'] > 0.75]
    for ppi in high_score.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1]
        afmult_score = high_score.loc[ppi, 'iptm+ptm']
        found = False
        for dir in afmult_dirs:
            prot1, prot2 = os.path.dirname(dir).split('/')[-1].split('__')[0], os.path.dirname(dir).split('/')[-1].split('__')[1]
            protA_dir, protB_dir = min([prot1, prot2]), max([prot1, prot2])
            if f"{protA}__{protB}" == f"{protA_dir}__{protB_dir}":
                os.system(f"cp {dir}ranked_0.pdb {filtered_files_out}{protA}__{protB}-{afmult_score*100:.0f}-ranked_0.pdb")


if False:
### Move high scoring HAF models
# vaccinia
    vacc_haf_het_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/heterodimer/"
    vacc_haf_ho_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/homodimer/"

    vacc_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
    vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial10.csv", index_col=0)
    high_score_vacc_HAF_out = f"{vacc_analysis_out}high_score_vacc_HAF_models/"
    os.makedirs(high_score_vacc_HAF_out, exist_ok=True)
    high_score_vac_HAF = vacc_master[vacc_master['Humphreys_maxcontact'] >= 0.8]
    for ppi in high_score_vac_HAF.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1]
        ppi_r = f"{protB}__{protA}"
        score = high_score_vac_HAF.loc[ppi, 'Humphreys_maxcontact']
        if os.path.exists(f"{vacc_haf_het_out}{ppi}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {vacc_haf_het_out}{ppi}_unrelaxed_model_3_ptm.pdb \
                      {high_score_vacc_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{vacc_haf_het_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {vacc_haf_het_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_vacc_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{vacc_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {vacc_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_vacc_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{vacc_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {vacc_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_vacc_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        else:
            print(f"{ppi} not found. skipping...")

# ASFV 
    asfv_haf_het_out = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/heterodimer/"
    asfv_haf_ho_out = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/homodimer/"

    asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
    asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial9.csv", index_col=0)
    high_score_asfv_HAF_out = f"{asfv_analysis_out}high_score_asfv_HAF_models/"
    os.makedirs(high_score_asfv_HAF_out, exist_ok=True)

    high_score_asfv_HAF = asfv_master[asfv_master['Humphreys_maxcontact'] >= 0.8]
    for ppi in high_score_asfv_HAF.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1]
        ppi_r = f"{protB}__{protA}"
        score = high_score_asfv_HAF.loc[ppi, 'Humphreys_maxcontact']
        if os.path.exists(f"{asfv_haf_het_out}{ppi}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {asfv_haf_het_out}{ppi}_unrelaxed_model_3_ptm.pdb \
                      {high_score_asfv_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{asfv_haf_het_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {asfv_haf_het_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_asfv_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{asfv_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {asfv_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_asfv_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        elif os.path.exists(f"{asfv_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb"):
            os.system(f'cp {asfv_haf_ho_out}{ppi_r}_unrelaxed_model_3_ptm.pdb \
                      {high_score_asfv_HAF_out}{ppi}_unrelaxed_model_3_ptm.pdb')
        else:
            print(f"{ppi} not found. skipping...")


