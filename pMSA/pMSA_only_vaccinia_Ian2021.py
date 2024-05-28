
import os, glob, pdb
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa


project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/Ian2021/"
combined_msa_out = f"{project_dir}combined_msas/"
pmsa_out = f"{project_dir}pmsas/"
temp = f"{project_dir}tmp/"
os.makedirs(pmsa_out, exist_ok=True)
#pair input from AFmult run
mccraith_afmult_fastas = glob.glob(f"/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/AFmult/*.fasta")

hhfilter = "/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter"
pairs = []
for fasta in mccraith_afmult_fastas: 
    descs, seqs = bpar.read_fasta(fasta)
    pairs.append((descs[0].split('__')[0], descs[1].split('__')[0]))
breakpoint()

if True:
    for pair in pairs:
        proteinA_msa = f"{combined_msa_out}{pair[0]}.fasta"
        proteinB_msa = f"{combined_msa_out}{pair[1]}.fasta"
        # generate paired and unpaired pMSAs and ouptut to .a3m files. This returns the raw data for downstream stats
        pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seq =  pmsa.generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99, unpaired_tr=90, temp=temp,  hhfilter=hhfilter)
