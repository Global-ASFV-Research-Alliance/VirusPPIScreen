import os, sys, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random

def generate_positive_controls(pairs, results, positive_protein_pairs):
    pos_ctrl = {}
    not_found_protein_pairs = set(positive_protein_pairs)
    for i in range(len(pairs)):
        for proteinA, proteinB in positive_protein_pairs:
            if proteinA in pairs[i] and proteinB in pairs[i]:
                pos_ctrl[pairs[i]] = results[i]
                if (proteinA, proteinB) in not_found_protein_pairs:
                    not_found_protein_pairs.remove((proteinA, proteinB))
                continue
    if not_found_protein_pairs:
        print(f"Not all positive controls found")
        for entry in not_found_protein_pairs:
            print(f"{entry[0]} {entry[1]} pair not found")
    else: 
        print(f"all positive control pairs found in results")
    no_value = 0
    delete = []
    for key in pos_ctrl:
        if pd.isna(pos_ctrl[key]):
            delete.append(key)
            no_value += 1
    for key in delete:
        del pos_ctrl[key]
    print(f"{100*no_value/len(pos_ctrl):.0f}% of control list have value of None\n{len(pos_ctrl)} total positive controls")
    return pos_ctrl, not_found_protein_pairs

def generate_negative_controls(pairs, results, pos_ctrl, ratio):
    """
    pairs is list of pMSA pair names, results is list of experimental values, pos_ctrl is dict from generate_positive_controls
    ratio is 1000 if there are 1:1000 true PPIs:total possible PPIs
    """
    fil_pairs, fil_results = [], []
    neg_ctrl = {}
    for i in range(len(pairs)):
            if pairs[i] not in pos_ctrl:
                fil_pairs.append(pairs[i])
                fil_results.append(results[i])
    # check to see if ratio doesn't work for number of pos controls
    num_neg_controls = int(ratio*len(pos_ctrl))
    total_pairs = len(pairs)
    if num_neg_controls > 0.8 * total_pairs:
        print(f"there are too many negative controls. {len(results)} pairs total. {num_neg_controls} requested for {len(pos_ctrl)} positive controls")
        downsize_ratio =  (0.8 * total_pairs) / num_neg_controls
        num_neg_controls = round(0.8 * total_pairs) # downsize to 80% of all data
        num_pos_controls = round(downsize_ratio * len(pos_ctrl))
        neg_indicies = random.sample(range(len(fil_pairs)), num_neg_controls)
        pos_indicies = random.sample(range(len(pos_ctrl)), num_pos_controls)
        for index in neg_indicies:
            neg_ctrl[fil_pairs[index]] = fil_results[index]
        pos_pairs, pos_results = [], []
        for key in pos_ctrl:
            pos_pairs.append(key)
            pos_results.append(pos_ctrl[key])
        pos_ctrl.clear()
        for index in pos_indicies:
            pos_ctrl[pos_pairs[index]] = pos_results[index]
        print(f'downsized to {num_pos_controls} positive controls and {num_neg_controls}') 
    else:                                     
        #filter lists of positive controls
       
        neg_ctrl = {}
        
        random_indices = random.sample(range(len(fil_pairs)), num_neg_controls)
        for index in random_indices:
            neg_ctrl[fil_pairs[index]] = fil_results[index]
    return pos_ctrl, neg_ctrl

def calculate_precision_recall(pos_ctrl, neg_ctrl, threshold):
    pos_above_tr = 0
    neg_above_tr = 0
    for key in pos_ctrl:
        if pos_ctrl[key] > threshold:
            pos_above_tr += 1
    for key in neg_ctrl:
        if neg_ctrl[key] > threshold:
            neg_above_tr += 1
    total_above_tr = pos_above_tr + neg_above_tr
    precision = pos_above_tr / total_above_tr
    recall = pos_above_tr / len(pos_ctrl)
    return precision, recall

def plot_line(X, Y, xlabel, ylabel, basename):
    # Define the font properties
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
           }

    plt.figure(figsize=(5, 5))  # Set the figure size
    plt.plot(X, Y, marker='o')  # Plot the line graph with circle markers

    # Set labels with the defined font properties
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)

    # Set axis limits
    plt.xlim(0,1)
    plt.ylim(0,1)

    # Set tick label sizes
    plt.xticks(fontsize=font['size'])
    plt.yticks(fontsize=font['size'])

    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{basename}_plt.png")
    plt.close()

def plot_scatter(X, Y, x_label, y_label, basename, plottype):
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
           }
    plt.figure(figsize=(5, 5))  # Set the figure size
    plt.scatter(X, Y, s=2)  # Create a scatter plot
    plt.xlim(0,10)
    plt.ylim(0,10)
    plt.xlabel(x_label, fontdict=font) #Set the x-axis label
    plt.ylabel(y_label, fontdict=font) # Set the y-axis label
    # Set tick label sizes
    plt.xticks(fontsize=font['size'])
    plt.yticks(fontsize=font['size'])
    plt.tight_layout()  # Adjust the layout to prevent clipping
    plt.savefig(f"{basename}_{plottype}.png")  # Save the plot to the specified file
    plt.close()  # Close the pl
#results_file, positive_db_file = sys.arg[1], sys.argv[2]
# input
results_file = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/output/jackhmmer_all_virus/ASFV_ASFV_clust99/ASFV_ASFV_clust99_PPIscores_combined_Neff90.csv"
positive_db_file_all = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_PPI_gold_and_polII.csv"
positive_db_gold = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_PPI_gold_stds.csv"
positive_db_polII = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_RNApolII_pairs.csv"
results_col = "RF2t++ scores"
ratio = 500 # ratio 
output_basename = 'output/test_recall'
# data prep
positive_db_df = pd.read_csv(positive_db_gold)
positive_protein_pairs = [(row['proteinA'], row['proteinB']) for _, row in positive_db_df.iterrows()]
results_df = pd.read_csv(results_file, index_col=0)
pairs = results_df.index.tolist()
results = results_df[results_col].tolist()
# generate controls 
pos_ctrl, not_found_protein_pairs = generate_positive_controls(pairs, results, positive_protein_pairs)
pos_ctrl, neg_ctrl = generate_negative_controls(pairs, results, pos_ctrl, ratio)
# generate precision recall curve
precision, recall = [], []
threshold_arr = np.arange(0, 1, .01)
for threshold in threshold_arr:
    p, r = calculate_precision_recall(pos_ctrl, neg_ctrl, threshold)
    precision.append(p)
    recall.append(r)
plot_line(recall, precision, 'Recall', 'Precision', output_basename)

#plot signal vs various parameters
RF2 = results_df['RF2t+ max contact probability'].tolist()
RF2plus = results_df['RF2t++ scores'].tolist()
Neff = results_df['Neff_90'].tolist()
Naln = results_df['#_alns'].tolist()
Nclust = results_df['Nseq_90'].tolist()
plot_scatter(Neff, RF2plus, 'Neff_90', 'RF2+ max contact prob', 'ASFV_clust99_90', 'Neff_vs_RF2')
plot_scatter(Naln, RF2, '# aln', 'RF2+ max contact prob', 'ASFV_clust99_90', 'NumAln_vs_RF2')
plot_scatter(Nclust, RF2, '# clust 90', 'RF2+ max contact prob', 'ASFV_clust99_90', 'Nclust80_vs_RF2')




