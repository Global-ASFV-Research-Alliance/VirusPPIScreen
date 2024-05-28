# this script is to generate the AlphaFold multimer input for the high scoring ASFV-ASFV heterodimer input
# that were found via Ian et al 2021 screening of the mega pMSA algorithm

import os, sys, glob, pdb
import pandas as pd
import common.bio_parsers as bpar

fasta_out = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/afmult_ian-tophits/"
vacc_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
descsi, seqsi = bpar.read_fasta(vacc_proteome)
descs, seqs = [], []

# clean up vaccinia naming
for desc, seq in zip(descsi, seqsi):
    name = desc.split('|')[1]
    descs.append(name)
    if seq[-1] == '*':
        seqs.append(seq[:-1])
    else:
        seqs.append(seq)

vacc_ian_heterodimer_dataset = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/vaccinia-WR_ian_mega_heterodimer.csv"

vacc_maxcontact = pd.read_csv(vacc_ian_heterodimer_dataset, index_col=0)

high_score = vacc_maxcontact[vacc_maxcontact['max_contact_prob'] >= 0.8] 
high_score_ppi = high_score.index.tolist()

if False:
    for ppi in high_score_ppi:
        proteinA, proteinB = ppi.split('__')[0], ppi.split('__')[1]
        seqA, seqB = '', ''
        for desc, seq in zip(descs, seqs):
            if proteinA == desc:
                seqA = seq
            if proteinB == desc:
                seqB = seq
        if seqA != '' and seqB != '':
            with open(f"{fasta_out}{ppi}.fasta", 'w') as f:
                f.write(f">{proteinA}\n{seqA}\n>{proteinB}\n{seqB}")

# split into 6 tranches
project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/afmult_ian-tophits/"

tranche1 = f"{project_dir}tranche1/"
tranche2 = f"{project_dir}tranche2/"
tranche3 = f"{project_dir}tranche3/"
tranche4 = f"{project_dir}tranche4/"
tranche5 = f"{project_dir}tranche5/"
tranche6 = f"{project_dir}tranche6/"

os.makedirs(tranche1, exist_ok=True), os.makedirs(tranche2, exist_ok=True), os.makedirs(tranche3, exist_ok=True)
os.makedirs(tranche4, exist_ok=True), os.makedirs(tranche5, exist_ok=True), os.makedirs(tranche6, exist_ok=True)

pmsas = glob.glob(f"{project_dir}*.fasta")
num_pmsas = len(pmsas)
print(f"there are {num_pmsas} pMSAs. Splitting into six tranches")
breakpoint()

 
if True:
    for i in range(0, int(num_pmsas/6)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche1}{basename}")
    for i in range(int(num_pmsas/6), int(num_pmsas/6 * 2)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche2}{basename}")
    for i in range(int(num_pmsas/6 * 2), int(num_pmsas/6 * 3)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche3}{basename}")
    for i in range(int(num_pmsas/6 * 3), int(num_pmsas/6 * 4)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche4}{basename}")
    for i in range(int(num_pmsas/6 * 4), int(num_pmsas/6 * 5)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche5}{basename}")
    for i in range(int(num_pmsas/6 * 5),int( num_pmsas)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche6}{basename}")

