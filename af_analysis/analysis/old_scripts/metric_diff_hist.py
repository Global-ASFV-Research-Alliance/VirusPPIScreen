import os, sys, math, statistics, glob, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict

def histogram(data, x_label, x_lim, outputdir, exp_tag, fontsize=16, annotation=(0.6, 0.6)):
    # input data as list of values
    # have RF2t++ x axis lim go to 2
    data = np.array(data) # convert to numpy array
    # Calculate statistics
    mean_value, median_value, stdev, num_values = np.nanmean(data), np.nanmedian(data), np.nanstd(data), data.size
    if num_values < 100:
        n_bins = 3 * round(np.sqrt(num_values))
    else: 
        n_bins = round(np.sqrt(num_values))
    plt.figure(figsize=(6, 6))
    plt.hist(data, bins=n_bins, edgecolor='black', alpha=0.7, density=False)
    plt.xlim(x_lim[0], x_lim[1])

    # Set font sizes using the fontsize variable
    plt.rcParams.update({'font.size': fontsize})
    plt.annotate(f"Mean: {mean_value:.1f}\nMedian: {median_value:.1f}\nStd Dev: {stdev:.1f}\nn: {num_values}\nn_bins: {n_bins}", xy=annotation, xycoords='axes fraction', fontsize=fontsize-6)
    plt.title(f"{exp_tag}_{x_label}", fontsize=fontsize)
    plt.xlabel(x_label+" diff", fontsize=fontsize)
    plt.ylabel('Counts', fontsize=fontsize)
    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)

    # Save the plot
    plt.tight_layout()
    plt.savefig(f"{outputdir}/{exp_tag}_{x_label.split(' ')[0]}diff.png")
    plt.close()

def generate_positive_controls(ctrl_db_file, results_file, results_col):
    positive_db_df = pd.read_csv(ctrl_db_file)
    positive_protein_pairs = [(row['proteinA'], row['proteinB']) for _, row in positive_db_df.iterrows()]
    results_df = pd.read_csv(results_file, index_col=0)
    pairs = results_df.index.tolist()
    results = results_df[results_col].tolist()
    pos_ctrl = {}
    not_found_protein_pairs = set(positive_protein_pairs)
    for i in range(len(pairs)):
        for proteinA, proteinB in positive_protein_pairs:
            if proteinA in pairs[i] and proteinB in pairs[i]:
                pos_ctrl[pairs[i].split('-pMSA')[0]] = results[i]
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
    print(f"{100*no_value/len(pos_ctrl):.0f}% of control list have value of None\n{len(pos_ctrl)} total positive controls")
    return pos_ctrl, not_found_protein_pairs

"""
usage: python3 scripts/metric_diff_hist.py /path/to/datafile1 /path/to/datafile2 /path/to/outputdir 'exp tag' 'column of interest'
python3 scripts/metric_diff_hist.py /Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_100/stats/jhmmer_ASFV-ASFV_Clust100_PPIscores_combined.csv /Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/output/jackhmmer_all_virus/ASFV_ASFV_clust99/ASFV_ASFV_clust99_PPIscores_combined.csv output/jackhmmer_all_virus 'ASFV_clust100vs99r2' 'RF2t+ max contact probability' 
"""

#data1_file, data2_file, output_dir, exp_tag, column_name = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
data1_file = '/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_100/stats/jhmmer_ASFV-ASFV_Clust100_PPIscores_combined.csv'
data2_file = '/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/output/jackhmmer_all_virus/ASFV_ASFV_clust99/ASFV_ASFV_clust99_PPIscores_combined.csv'
output_dir, exp_tag, column_name = 'output/jackhmmer_all_virus', 'ASFV_clust100vs99r2', 'RF2t+ max contact probability'
df1, df2 = pd.read_csv(data1_file, index_col=0), pd.read_csv(data2_file, index_col=0)
diff, omit = [], 0
"""
below will only work for -pMSAfilexact vs -pMSAfil
"""
for index in df1.index:
    index2_a = f"{index.split('-pMSA')[0]}-pMSAfil"
    index2_b = f"{index.split('-pMSA')[0].split('__')[1]}__{index.split('-pMSA')[0].split('__')[0]}-pMSAfil"
    try:
        diff.append(df1.loc[index, column_name]-df2.loc[index2_a, column_name])
    except KeyError:
        try:
            diff.append(df1.loc[index, column_name]-df2.loc[index2_b, column_name])
        except KeyError:
            omit += 1
            continue       
    except:
        pdb.set_trace()
print(f"There are {omit} keys omitted")

x_lim = (min(diff), max(diff))
x_label = 'RF2t+'
#pdb.set_trace()
histogram(diff, x_label, x_lim, output_dir, exp_tag, fontsize=20)

#control section
positive_db_file_all = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_PPI_gold_and_polII.csv"
positive_db_gold = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_PPI_gold_stds.csv"
positive_db_polII = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/data/control_sets/ASFV_RNApolII_pairs.csv"

pos_ctrl1, not_found_protein_pairs1 = generate_positive_controls(positive_db_gold, data1_file, column_name)
pos_ctrl2, not_found_protein_pairs2 = generate_positive_controls(positive_db_gold, data2_file, column_name)
ctrl_diff = []
#pdb.set_trace()
for key in pos_ctrl1:
    try:
        ctrl_diff.append(pos_ctrl1[key]-pos_ctrl2[key])
    except KeyError:
        alt_key = f"{key.split('__')[1]}__{key.split('__')[0]}"
        ctrl_diff.append(pos_ctrl1[key]-pos_ctrl2[alt_key])
x_lim = (min(ctrl_diff), max(ctrl_diff))
histogram(ctrl_diff, x_label, x_lim, output_dir, f"{exp_tag}_pos-ctrl_gold", fontsize=20, annotation=(.2, 0.6))
