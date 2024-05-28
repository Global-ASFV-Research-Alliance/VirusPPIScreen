# These defs process RoseTTAFold 2-track results according to Humphrey's et al 2021 methods
import os, sys, glob, pdb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def read_NPZ_to_arr(npz_path):
# this def inputs a single path to a npz file, checks if there is only one array contained 
# and returns the numpy array. This is specific to the RoseTTAfold 2-track, RF2t, .npz output 
    npz = np.load(npz_path)
    if len(npz.keys()) == 1: #make sure npz input has only 1 array
        arr = npz[npz.files[0]]
    else:
        print('Error: Input .npz has more than one array') 
    return arr

def RF2t_plus_termini_correction(arr):
# this definition inputs a RF2t numpy array of contact probabilities and excludes the last ten amino acids
# of the first protein and the first ten amino acids of the second protein when calculating the max contact PPI
    arr[-10:, :] = 0 # set the last 10 residues on C terminus of protein 1 to zero
    arr[:, :10] = 0 # set the first 10 residues on N terminus of protein 2 to zero
    max_contact = np.max(arr)
    return max_contact

def APC_RF2t_plusplus_score(ppis, max_contact_probs, sym=True):
# this is implimenting an average product correction, APC, to the max_contact probability scores
# genes1 and genes2 are the list of genes in pMSAs. They are identical for ASFV-ASFV pMSAs
# lists file_basenames and max_contact_probs must have the same indicies
    if '__' in ppis[0]: #switch to process genes with the format gene1-AA#__gene2-AA#
        # get list of proteins from file names
        if sym: #this extracts the protein list assuming symettric pairing, ie ASFV__ASFV pMSAs (all vs all)
            proteins, pairs = [], []
            for pair in ppis:
                X, Y, = pair.split('__')[0], pair.split('__')[1] #extract the two gene names
                if (X,Y) in pairs or (Y,X) in pairs: #make sure you don't accidently load duplicates 
                    print(f"Input .npz contain duplicate pairs. ex. {pair}.npz is present more than once. APC not valid, exiting")
                    exit()
                else: 
                    pairs.append((X,Y))
                if X not in proteins:
                    proteins.append(X)
                if Y not in proteins:
                    proteins.append(Y)
            score_matrix = pd.DataFrame(data=float(0), index=proteins, columns=proteins) #initialize df with zeros
            # Generate Prob X,Y matrix
            for pair, prob in zip(ppis, max_contact_probs):
                X, Y, = pair.split('__')[0], pair.split('__')[1] #extract the two gene names
                score_matrix.loc[X, Y] = prob
                score_matrix.loc[Y, X] = prob
        else: #this extracts the protein1 and protein2 list if asymetric, ie. ASFV__Sscrofa pMSAs
            proteins1, proteins2, pairs = [], [], []
            for pair in ppis:
                X, Y, = pair.split('__')[0], pair.split('__')[1] #extract the two gene names
                if (X,Y) in pairs: #make sure you don't accidently load duplicates 
                    print(f"Input .npz contain duplicate pairs. ex. {pair}.npz is present more than once. APC not valid, exiting")
                    exit()
                else: 
                    pairs.append((X,Y))
                if X not in proteins1:
                    proteins1.append(X)
                if Y not in proteins2:
                    proteins2.append(Y)
            score_matrix = pd.DataFrame(data=float(0), index=proteins1, columns=proteins2) #initialize df with zeros
            # Generate Prob X,Y matrix
            for pair, prob in zip(ppis, max_contact_probs):
                X, Y, = pair.split('__')[0], pair.split('__')[1] #extract the two gene names
                score_matrix.loc[X, Y] = prob
        # calculate series of X, Y vs all, and all vs all value
        protein1_vs_all = score_matrix.sum(axis=1) 
        protein2_vs_all = score_matrix.sum(axis=0)
        all_vs_all = protein1_vs_all.sum()
        scores = [] #APC score of each PPI probability
        for pair, prob in zip(ppis, max_contact_probs):
            X, Y, = pair.split('__')[0], pair.split('__')[1] #extract the two gene names
            if prob >= 0.9:
                score = prob + prob - protein1_vs_all[X] * protein2_vs_all[Y] / all_vs_all
            elif prob < 0.9:
                score = prob - protein1_vs_all[X] * protein2_vs_all[Y] / all_vs_all
            scores.append(score)
    else:
        print('Incompatible .npz format, script exited')
        exit() 
    return scores

def process_APC_RF2t_plusplus_npz (npz_glob, sym=True, name_switch=None):
# high level def that processes .npz files from RF2t model. Conducts N&C terimini corrections (RF2t+) and APC scoring (RF2t++)
# input file path glob pattern string, npz_glob, the output directory string, outputdir, the experiment_name string for summary file naming, 
# and if the set of pMSAs to generate the npz files are symettric or not. sym=True would be ASFV__ASFV all vs all pMSAs where 
# sym=False would be ASFV_Sscrufa pMSAs. Default sym is True.
    ppis, max_contact_probs_RF2t_plus = [], [] #lists to hold the PPI names and corresponding max contact probabilities
    for npz_path in glob.iglob(npz_glob):
        basename = os.path.basename(npz_path).split('.')[0]
        chain1 = int(basename.split('__')[0].split('-AA')[1])
        chain2 = int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if name_switch == None:
            gene1 = basename.split('-AA')[0] 
            gene2 = basename.split('__')[1].split('-AA')[0]
        elif name_switch == 'fmdv':
            gene1 = basename.split('__')[0].split('768_')[1].split('-AA')[0]
            gene2 = basename.split('__')[1].split('768_')[1].split('-AA')[0]
        protA, protB = min([gene1, gene2]), max([gene1, gene2])
        ppi = f"{protA}__{protB}"
        arr = read_NPZ_to_arr(npz_path)
        max_RF2t_plus = RF2t_plus_termini_correction(arr)
        ppis.append(ppi)
        max_contact_probs_RF2t_plus.append(max_RF2t_plus)
    scores_RF2t_plusplus = APC_RF2t_plusplus_score(ppis, max_contact_probs_RF2t_plus, sym)
    # Export all max contact probabilites to .csv
    df = pd.DataFrame(data={'RF2t+ max contact probability':max_contact_probs_RF2t_plus, 'RF2t++ scores':scores_RF2t_plusplus},
                      index=ppis)
    return df

def process_RF2t_plus_npz(npz_glob, name_switch=None):
# high level def that processes .npz files from RF2t model. Conducts N&C terimini corrections (RF2t+) and APC scoring (RF2t++)
# input file path glob pattern string, npz_glob, the output directory string, outputdir, the experiment_name string for summary file naming, 
# and if the set of pMSAs to generate the npz files are symettric or not. sym=True would be ASFV__ASFV all vs all pMSAs where 
# sym=False would be ASFV_Sscrufa pMSAs. Default sym is True.
    ppis, max_contact_probs_RF2t_plus = [], [] #lists to hold the PPI names and corresponding max contact probabilities
    for npz_path in glob.iglob(npz_glob):
        basename = os.path.basename(npz_path).split('.')[0]
        chain1 = int(basename.split('__')[0].split('-AA')[1])
        chain2 = int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if name_switch == None:
            gene1 = basename.split('-AA')[0] 
            gene2 = basename.split('__')[1].split('-AA')[0]
        elif name_switch == 'fmdv':
            gene1 = basename.split('__')[0].split('768_')[1].split('-AA')[0]
            gene2 = basename.split('__')[1].split('768_')[1].split('-AA')[0]
        protA, protB = min([gene1, gene2]), max([gene1, gene2])
        ppi = f"{protA}__{protB}"
        arr = read_NPZ_to_arr(npz_path)
        max_RF2t_plus = RF2t_plus_termini_correction(arr)
        ppis.append(ppi)
        max_contact_probs_RF2t_plus.append(max_RF2t_plus)
    # Export all max contact probabilites to .csv
    df = pd.DataFrame(data={'RF2t+ max contact probability':max_contact_probs_RF2t_plus},
                      index=ppis)
    return df


def rf2t_npz_analysis(npz_file_dir, symmetry=True, APC=True, name_switch=None):
# Usage: $ python3 RF2t_npz_analysis.py <input *.npz file directory> <output directory> <experiment name> symmetry=<True/False> APC=<True/False>
# Usage Example (order matters, must input all): $ python3 RF2t_npz_analysis.py /path/to/data /path/to/output test_experiment symmetry=True APC=True
# python3 scripts/RF2t_npz_analysis.py /Users/jacobfenster/Documents/ASFV_Postdoc/Data/Internal/RoseTTAfold_results/jackhmmer_all_virus/ASFV_ASFV_pMSA_clust100_Output output/jackhmmer_all_virus/ASFV_ASFV_clust100 jhmmer_ASFV-ASFV_Clust100 symmetry=True APC=True
    npz_glob = npz_file_dir + '*.npz'
    if APC:
        df = process_APC_RF2t_plusplus_npz(npz_glob, symmetry, name_switch=name_switch)
    else:
        df = process_RF2t_plus_npz(npz_glob, name_switch=name_switch)
    print(f'RF2t_npz_analysis complete for dir {npz_file_dir}')
    return df
