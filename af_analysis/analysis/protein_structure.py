# This contains functions to impliment protein strucutre analysis and comparason of model outputs
import os, sys, subprocess, pdb
import pandas as pd
import numpy as np

lddt_container = "/lustrefs/software/public/containers/openstructure-2.7.0-JF.sif"

reference = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/high_score_heterodimer/EP296R__DP96R/ranked_0.pdb"
model = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/AFIan2021/asfv-asfv_homo_heterodimers/heterodimer/EP296R__DP96R_unrelaxed_model_3_ptm.pdb" 

# implimentation of lDDT
def lddt(reference, model):
    """
    Singularity implementation of lDDT calculation
    singularity run --app lDDT <IMAGE> [options] <mod1> [mod1 [mod2]] <re1>[,ref2,ref3]

    """
    command = f"singularity run --app lDDT {lddt_container} {model} {reference}"
    breakpoint()
    result = subprocess.run(f"{command} > lDDT.out", shell=True, text=True, capture_output=True)
    if result.returncode == 0:
        x=1

def ost(reference, model):
    """
    runs the OpenStructure OST in singularity container
    """
    command = f"singularity run --app OST {lddt_container} -m {model} -r {reference}"
    breakpoint()
    result = subprocess.run(f"{command}", shell=True, text=True, capture_output=True)

# implementation of TM

ost(reference, model)
