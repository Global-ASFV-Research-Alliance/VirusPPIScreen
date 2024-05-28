# this script is to generate the AlphaFold multimer input for the high scoring ASFV-ASFV heterodimer input
# that were found via Ian et al 2021 screening of the mega pMSA algorithm

import os, sys, glob, pdb
import pandas as pd
import common.bio_parsers as bpar

fasta_out = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/alphafold_multimer/mid_high_score_heterodimer/"
asfv_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/Georgia-2007-LR743116_trans-rename.fa"
descs, seqs = bpar.read_fasta(asfv_proteome)

asfv_heterodimer_Ian2021_dataset = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/asfv-asfv_heterodimer_maxcontact.csv"

asfv_maxcontact = pd.read_csv(asfv_heterodimer_Ian2021_dataset, index_col=0)

high_score = asfv_maxcontact[asfv_maxcontact['max_contact_prob'] >= 0.9] 
mid_high_score = asfv_maxcontact[(asfv_maxcontact['max_contact_prob'] < 0.9) & (asfv_maxcontact['max_contact_prob'] >= 0.8)] 
# altered to look at those between 0.8 and 0.9 on 3/1/24
breakpoint()

high_score_ppi = mid_high_score.index.tolist()

for ppi in high_score_ppi:
    proteinA, proteinB = ppi.split('__')[0], ppi.split('__')[1]
    seqA, seqB = '', ''
    for desc, seq in zip(descs, seqs):
        if proteinA == desc:
            seqA = seq[:-1]
        if proteinB == desc:
            seqB = seq[:-1]
    if seqA != '' and seqB != '':
        with open(f"{fasta_out}{ppi}.fasta", 'w') as f:
            f.write(f">{proteinA}\n{seqA}\n>{proteinB}\n{seqB}")

