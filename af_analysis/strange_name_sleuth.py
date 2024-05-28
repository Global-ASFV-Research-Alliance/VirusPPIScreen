# This script is to figure out what is going on with the strainge names in the mega pMSAs
import os, sys, glob, json, pdb
import pickle
import numpy as np
import pandas as pd 
# new to the AFppi 
import subprocess, math, time
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar


def _strange_name_sleuth(pmsa):
    basename = os.path.basename(pmsa).split('.')[0]
    gene1 = basename.split('-AA')[0] 
    gene2 = basename.split('__')[1].split('-AA')[0]
    chain1 = int(basename.split('__')[0].split('-AA')[1])
    chain2 = int(basename.split('__')[1].split('-AA')[1].split('-')[0])
    descs, seqs = bpar.read_fasta(pmsa)
    paired_descs, paired_seqs = [], []
    chain1_descs, chain1_seqs = [], []
    chain2_descs, chain2_seqs = [], []
    for desc, seq in zip(descs, seqs):
        if '__' in desc and len(desc.split('__')[0].split('.')) == 2 and len(desc.split('__')[1].split('.')) == 2:
            paired_descs.append(desc), paired_seqs.append(seq)
        elif all(char == '-' for char in seq[chain1:]):
            chain1_descs.append(desc), chain1_seqs.append(seq[:chain1])
        elif all(char == '-' for char in seq[:chain1]):
            chain2_descs.append(desc), chain2_seqs.append(seq[chain1:])
        else:
            breakpoint()
            x=1

fmdv_het_pmsa_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/AFIan2021/mega_pMSA/pmsas/paired_unpaired/"
fmdv_het_pmsas = glob.glob(f"{fmdv_het_pmsa_dir}*.a3m")
for pmsa in fmdv_het_pmsas:
    _strange_name_sleuth(pmsa)

