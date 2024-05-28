import os, sys, glob, pdb
import pandas as pd
import numpy as np
import pickle
import analysis.AFppi as ppi
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar
sys.path.append("/home/jacob.fenster/scripts/alphafold-2.3.2/")
from alphafold.data import parsers
from alphafold.data import templates
from alphafold.data import mmcif_parsing
from Bio.PDB import MMCIF2Dict
import time


# This script is to conduct hhsearch of the pdb70 and analyze templates that will go into the AF pipeline
# Looks at the pdb release date to see if the structure would have been used during training
# this script is to be run inside the singularity container /lustrefs/software/public/containers/alphafold_2.3.2-1.sif
# AF multimer output dir
mccraith_afmult_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/AFmult/"
mccraith_afmult_fastas = glob.glob(f"/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/*.fasta")
asfv_gold_afmult = "/lustrefs/fadru/projects/asfv-ppi/Output/20231204_AFtests/ASFV_ASFV-Gold/" 

# software input 
asfv_output = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/ASFV/"
vaccinia_output = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/Vaccinia/"
temp = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/tmp/"
hhsearch = "hhsearch"
reformat = "/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl"
# parse pdb release dates from in house created .txt file
pdb_release = templates._parse_release_dates("/home/jacob.fenster/scripts/af_analysis/data/pdb_release_date.txt")

mccraith_pairs = [] # holds McCraith PPI pairs parsed from AF mult run
mccraith_proteins = set()
# extract PPI pairs 
for fasta in mccraith_afmult_fastas: 
    descs, seqs = bpar.read_fasta(fasta)
    mccraith_pairs.append((descs[0].split('__')[0], descs[1].split('__')[0]))
    mccraith_proteins.add(descs[0].split('__')[0]), mccraith_proteins.add(descs[1].split('__')[0])
breakpoint()
# conduct template search for McCraith proteins
for protein_name in mccraith_proteins:
    hrr_string = ppi.afmult_template_search(protein_name, mccraith_afmult_out, temp=temp, reformat='/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl', hhsearch='hhsearch', pdb70='/lustrefs/public/databases/alphafold/pdb70/')
    template_hits = parsers.parse_hhr(hrr_string)
    # Error finding the pdb ID in the pdb_release data. Need to troubleshoot
    # also need to find a way to rank the hits

