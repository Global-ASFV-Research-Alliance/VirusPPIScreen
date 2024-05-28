import os as os
import glob
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import pdb

def histogram(data, x_label, x_lim, outputdir, experiment_name):
# input data as list of values
# have RF2t++ x axis lim go to 2
    data = np.array(data) #convert to numpy array
    # Calculate statistics
    mean_value, median_value, stdev, num_values = np.mean(data), np.median(data), np.std(data), data.size
    if num_values < 100:
        n_bins = 3 * int(np.sqrt(num_values))
    else: 
        n_bins = int(np.sqrt(num_values))
    plt.figure(figsize=(6, 6))
    plt.hist(data, bins=n_bins, edgecolor='black', alpha=0.7, density=False)
    plt.xlim(0, x_lim)
    plt.rcParams.update({'font.size': 20})
    plt.annotate(f"Mean: {mean_value:.1f}\nMedian: {median_value:.1f}\nStd Dev: {stdev:.1f}\nn: {num_values}\nn_bins: {n_bins}", xy=(0.6, 0.8), xycoords='axes fraction', fontsize=12)
    plt.title(experiment_name)
    plt.xlabel(x_label)
    plt.ylabel('Counts')
    # Save the plot
    plt.tight_layout()
    plt.savefig(outputdir+experiment_name+x_label+'.png')
    plt.close()
    return

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

def APC_RF2t_plusplus_score(file_basenames, max_contact_probs, sym=True):
# this is implimenting an average product correction, APC, to the max_contact probability scores
# genes1 and genes2 are the list of genes in pMSAs. They are identical for ASFV-ASFV pMSAs
# lists file_basenames and max_contact_probs must have the same indicies
    if '__' in file_basenames[0]: #switch to process genes with the format gene1-AA#__gene2-AA#
        # get list of proteins from file names
        if sym: #this extracts the protein list assuming symettric pairing, ie ASFV__ASFV pMSAs (all vs all)
            proteins, pairs = [], []
            for pair in file_basenames:
                X, Y, = pair.split('__')[0].split('-AA')[0], pair.split('__')[1].split('-AA')[0] #extract the two gene names
                if (X,Y) in pairs or (Y,X) in pairs: #make sure you don't accidently load duplicates 
                    print(f"Input .npz contain duplicate pairs. ex. {pair}.npz is present more than once. APC not valid, exiting")
                    exit()
                else: 
                    pairs.append((X,Y))
                if X not in proteins:
                    proteins.append(X)
                if Y not in proteins:
                    proteins.append(Y)
            score_matrix = pd.DataFrame(data=0, index=proteins, columns=proteins) #initialize df with zeros
            # Generate Prob X,Y matrix
            for pair, prob in zip(file_basenames, max_contact_probs):
                X, Y, = pair.split('__')[0].split('-AA')[0], pair.split('__')[1].split('-AA')[0] #extract the two gene names
                score_matrix.loc[X, Y] = prob
                score_matrix.loc[Y, X] = prob
        else: #this extracts the protein1 and protein2 list if asymetric, ie. ASFV__Sscrofa pMSAs
            proteins1, proteins2, pairs = [], [], []
            for pair in file_basenames:
                X, Y, = pair.split('__')[0].split('-AA')[0], pair.split('__')[1].split('-AA')[0] #extract the two gene names
                if (X,Y) in pairs: #make sure you don't accidently load duplicates 
                    print(f"Input .npz contain duplicate pairs. ex. {pair}.npz is present more than once. APC not valid, exiting")
                    exit()
                else: 
                    pairs.append((X,Y))
                if X not in proteins1:
                    proteins1.append(X)
                if Y not in proteins2:
                    proteins2.append(Y)
            score_matrix = pd.DataFrame(data=0, index=proteins1, columns=proteins2) #initialize df with zeros
            # Generate Prob X,Y matrix
            for pair, prob in zip(file_basenames, max_contact_probs):
                X, Y, = pair.split('__')[0].split('-AA')[0], pair.split('__')[1].split('-AA')[0] #extract the two gene names
                score_matrix.loc[X, Y] = prob
        # calculate series of X, Y vs all, and all vs all value
        protein1_vs_all = score_matrix.sum(axis=1) 
        protein2_vs_all = score_matrix.sum(axis=0)
        all_vs_all = protein1_vs_all.sum()
        scores = [] #APC score of each PPI probability
        for pair, prob in zip(file_basenames, max_contact_probs):
            X, Y, = pair.split('__')[0].split('-AA')[0], pair.split('__')[1].split('-AA')[0] #extract the two gene names
            if prob >= 0.9:
                score = prob + prob - protein1_vs_all[X] * protein2_vs_all[Y] / all_vs_all
            elif prob < 0.9:
                score = prob - protein1_vs_all[X] * protein2_vs_all[Y] / all_vs_all
            scores.append(score)
    else:
        print('Incompatible .npz format, script exited')
        exit() 
    return scores 

def process_APC_RF2t_plusplus_npz (npz_glob, outputdir, experiment_name, sym=True):
# high level def that processes .npz files from RF2t model. Conducts N&C terimini corrections (RF2t+) and APC scoring (RF2t++)
# input file path glob pattern string, npz_glob, the output directory string, outputdir, the experiment_name string for summary file naming, 
# and if the set of pMSAs to generate the npz files are symettric or not. sym=True would be ASFV__ASFV all vs all pMSAs where 
# sym=False would be ASFV_Sscrufa pMSAs. Default sym is True.
    file_basenames, max_contact_probs_RF2t_plus = [], [] #lists to hold the PPI names and corresponding max contact probabilities
    for npz_path in glob.iglob(npz_glob):
        basename = os.path.basename(npz_path).split('.')[0]
        arr = read_NPZ_to_arr(npz_path)
        max_RF2t_plus = RF2t_plus_termini_correction(arr)
        file_basenames.append(basename)
        max_contact_probs_RF2t_plus.append(max_RF2t_plus)
    scores_RF2t_plusplus = APC_RF2t_plusplus_score(file_basenames, max_contact_probs_RF2t_plus, sym)
    # Export all max contact probabilites to .csv
    df = pd.DataFrame(data={'PPI name':file_basenames, 'RF2t+ max contact probability':max_contact_probs_RF2t_plus, 'RF2t++ scores':scores_RF2t_plusplus})
    df.to_csv(outputdir+experiment_name+'_PPIscores.csv', index=0)
    histogram(max_contact_probs_RF2t_plus, 'RF2t+ max contact probability', 1, outputdir, experiment_name)
    histogram(scores_RF2t_plusplus, 'RF2t++ APC crrected scores', 2, outputdir, experiment_name)
    return

def process_RF2t_plus_npz(npz_glob, outputdir, experiment_name):
# high level def that processes .npz files from RF2t model. Conducts N&C terimini corrections (RF2t+) and APC scoring (RF2t++)
# input file path glob pattern string, npz_glob, the output directory string, outputdir, the experiment_name string for summary file naming, 
# and if the set of pMSAs to generate the npz files are symettric or not. sym=True would be ASFV__ASFV all vs all pMSAs where 
# sym=False would be ASFV_Sscrufa pMSAs. Default sym is True.
    file_basenames, max_contact_probs_RF2t_plus = [], [] #lists to hold the PPI names and corresponding max contact probabilities
    for npz_path in glob.iglob(npz_glob):
        basename = os.path.basename(npz_path).split('.')[0]
        arr = read_NPZ_to_arr(npz_path)
        max_RF2t_plus = RF2t_plus_termini_correction(arr)
        file_basenames.append(basename)
        max_contact_probs_RF2t_plus.append(max_RF2t_plus)
    # Export all max contact probabilites to .csv
    df = pd.DataFrame(data={'PPI name':file_basenames, 'RF2t+ max contact probability':max_contact_probs_RF2t_plus})
    df.to_csv(outputdir+experiment_name+'_PPIscores.csv', index=0)
    histogram(max_contact_probs_RF2t_plus, 'RF2t+ max contact probability', 1, outputdir, experiment_name)
    return


# Usage: $ python3 RF2t_npz_analysis.py <input *.npz file directory> <output directory> <experiment name> symmetry=<True/False> APC=<True/False>
# Usage Example (order matters, must input all): $ python3 RF2t_npz_analysis.py /path/to/data /path/to/output test_experiment symmetry=True APC=True
# python3 scripts/RF2t_npz_analysis.py /Users/jacobfenster/Documents/ASFV_Postdoc/Data/Internal/RoseTTAfold_results/jackhmmer_all_virus/ASFV_ASFV_pMSA_clust100_Output output/jackhmmer_all_virus/ASFV_ASFV_clust100 jhmmer_ASFV-ASFV_Clust100 symmetry=True APC=True
if len(sys.argv) != 6:
    print('Incorrect number of input arguments. Usage: $ python3 RF2t_npz_analysis.py <input *.npz file directory> <output directory> <experiment name> symmetry=<True/False> APC=<True/False>\nUsage Example: $ python3 RF2t_npz_analysis.py /path/to/data /path/to/output test_experiment symmetry=True APC=True')
    exit()
npz_file_dir, outputdir, experiment_name, symmetry, APC = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split('=')[1], sys.argv[5].split('=')[1]
if symmetry == 'True':
    symmetry = True
elif symmetry == 'False':
    symmetry = False
else:
    print('Invalid input for symmetry. Input symmetry=<True/False>. symmetric=True means all against all pMSAs as in ASFV__ASFV. symmetric=False would be ASFV__Sscrofa pMSAs')
    exit()
if APC != 'True' and APC != 'False':
    print('Invalid input for APC. Input APC=<True/False> to use or not use average product correction.')
    exit() 
npz_glob = npz_file_dir + '/*.npz'
outputdir = outputdir + '/'
if APC == 'True':
    process_APC_RF2t_plusplus_npz(npz_glob, outputdir, experiment_name, symmetry)
if APC == 'False':
    process_RF2t_plus_npz(npz_glob, outputdir, experiment_name)
print('RF2t_npz_analysis done')


