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
    plt.hist(data, bins=n_bins, edgecolor='black', alpha=0.7, density=True)
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

tax_data_file = "/Users/jacobfenster/Documents/ASFV_Postdoc/Entrez_API/output/jackhmmer_all_virus/ASFV_clust99/stats/pMSA_tax_summary.csv"
tax_data = pd.read_csv(tax_data_file, index_col=0)
per_ASFV = tax_data["%African swine fever virus"].tolist()
per_error_species = tax_data["%_error_species"].tolist()
per_ASFV_known = []
for pA, pe in zip(per_ASFV, per_error_species):
    per_ASFV_known.append(pA+pe)
x_label = '% ASFV of known'
x_lim = (0,100)
outputdir = '/Users/jacobfenster/Documents/ASFV_Postdoc/Entrez_API/output/jackhmmer_all_virus/ASFV_clust99/stats/figs'
exp_tag = 'pMSA_tax_density'
histogram(per_ASFV_known, x_label, x_lim, outputdir, exp_tag, fontsize=16, annotation=(0.6, 0.6))