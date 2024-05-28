import os, sys, math, statistics, glob, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict


def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of descs and a list of seqs.
    """
    descs = []
    seqs = []
    sequence = ""

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespaces
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):  # Description line
                if sequence:
                    seqs.append(sequence)
                    sequence = ""
                descs.append(line[1:])
            else:
                sequence += line

        if sequence:
            seqs.append(sequence)

    return descs, seqs

def query_ID_cov_gap_score(aln, query_aln):
    """ 
    input the alnignment and query alignment in fasta format and output %ID, %Coverage, and %gap score
    %ID = 100 * identical columns/aligned columns, %Coverage = 100 * aligned columns/query length
    %gap = 100 * gaps/query length where query length is number of gap-less columns in query (ignores insertion columns in target)
    """
    aligned, identical, gaps, query_len = 0, 0, 0, 0
    for a, b in zip(aln, query_aln):
        if a != '-' and b != '-':  # Skip gaps 
            aligned += 1
            query_len += 1
            if a == b:  # Count identical positions
                identical += 1
        elif a == '-' and b != '-': #count gaps when present in target column and not in query column
                gaps += 1
                query_len += 1 #query length does not count insertion rows in target
    if aligned == 0 and query_len != 0:
        return 0.0, 100*(aligned/query_len), 100*(gaps/query_len)
    elif aligned == 0 and query_len == 0:
        return 0.0, 0.0, 0.0 
    return 100*(identical/aligned), 100*(aligned/query_len), 100*(gaps/query_len)

def residue_coverage_values(alns, query_aln):
    """
    This one inputs a list of fasta formatted sequence alignments and outputs number of sequences that have residues
    at each position in the query, first aln
    """
    aln_len = len(query_aln)
    residues_cov_x = [i for i in range(aln_len)] 
    residues_cov_y = [0 for i in range(aln_len)]
    for aln in alns:
        for c in range(aln_len):
            if aln[c].isalpha():
                residues_cov_y[c] += 1
    return residues_cov_x, residues_cov_y

def combined_hist_scatter_plt(IDs_perct, covs_perct, residues_cov_x, residues_cov_y, Neff, basename):
    tick_fontsize, label_fontsize = 18, 20
    num_values = len(IDs_perct)
    if num_values < 100:
        n_bins = 3 * int(math.sqrt(num_values))
    else: 
        n_bins = int(math.sqrt(num_values))
    # Create a figure and a 1x2 grid of subplots
    fig, axs = plt.subplots(1, 2, figsize=(12, 6))
    # Add an overall figure title
    title_basename = basename.split('-pMSA')[0]
    fig.suptitle(f'{title_basename} Neff={Neff:.1f} Ntotal={num_values}', fontsize=label_fontsize)
    # Left subplot: Histogram
    axs[0].hist(IDs_perct, bins=n_bins, alpha=0.5, label='%ID', color='b', density=True)
    axs[0].hist(covs_perct, bins=n_bins, alpha=0.6, label='%Coverage', color='#F28C28', density=True)
    axs[0].set_xlabel('%ID or %Coverage', fontsize=label_fontsize)
    axs[0].set_ylabel('Frequency', fontsize=label_fontsize)
    axs[0].set_xlim([0, 100])
    axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), frameon=False, fontsize=(tick_fontsize-4), framealpha=0)
    # Right subplot: Scatter Plot
    axs[1].plot(residues_cov_x, residues_cov_y, 'k', linewidth=1, label='Residue Coverage')
    axs[1].set_xlabel('Position', fontsize=label_fontsize)
    axs[1].set_ylabel('Sequences', fontsize=label_fontsize)
    axs[1].set_xlim([0, max(residues_cov_x)])
    axs[1].set_ylim([0, max(residues_cov_y)+20])
    #update axis tick fonts
    for tick in axs[0].xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[0].yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[1].xaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    for tick in axs[1].yaxis.get_major_ticks():
        tick.label1.set_fontsize(tick_fontsize)
    # Remove top and right spines for left subplot
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    # Increase linewidth and tick size
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].tick_params(axis='both', width=2, size=6)
    axs[1].spines['left'].set_linewidth(2)
    axs[1].spines['bottom'].set_linewidth(2)
    axs[1].tick_params(axis='both', width=2, size=6)
    # Show the plot
    plt.rcParams.update({'font.size': 24})
    plt.tight_layout() #rect=[0, 0, 1, 1]
    plt.subplots_adjust(top=0.75)
    if workaround:
        plt.savefig(f"{output_dir}/{basename}_plt.png") #changed for workaround
    else:
        plt.savefig(f"{basename}_plt.png")
    plt.close()

def MSA_summary(MSA_list, Nseq_dict):
    """
    Inputs a list of paths to MSAs or pMSAs in fasta format, calculates query specific stats, Neff sequences, and adds to a row in a growing dataframe
    First line of MSA/pMSA must be query alignment.
    with the MSA speicifc stats it generates a histogram of the %ID and %Coverage distribution with summary stats.  
    """
    summary = pd.DataFrame(columns=['#_alns', 'Nseq_80', 'Neff_80', 'mean_ID', 'stdev_ID', 'mean_cov', 'stdev_cov'])
    for MSA in MSA_list:
        index = os.path.basename(MSA).split('.')[0] #name the summary row by the MSA/pMSA basename
        descs, alns = read_fasta(MSA)
        IDs, covs = [], []
        query_aln = alns[0] #first entry must be query sequence
        for i in range(len(alns)): #find vs query summary stats
            ID_perct, cov_perct, gap_perct = query_ID_cov_gap_score(alns[i], query_aln)
            IDs.append(ID_perct)
            covs.append(cov_perct)
        residues_cov_x, residues_cov_y = residue_coverage_values(alns, query_aln)
        if 'filexact' in index:
            Nseq_index = f"{index[:-8]}_fil80"
        else:
            Nseq_index = f"{index[:-3]}_fil80"
        Nseq, Neff = Nseq_dict[Nseq_index], Nseq_dict[Nseq_index]/math.sqrt(len(query_aln))
        summary.loc[index, '#_alns'], summary.loc[index, 'Nseq_80'], summary.loc[index, 'Neff_80'] = len(alns), Nseq, Neff
        summary.loc[index, 'mean_ID'], summary.loc[index, 'mean_cov'] = statistics.mean(IDs), statistics.mean(covs)
        if len(IDs) < 2:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = 0, 0
        else:
            summary.loc[index, 'stdev_ID'], summary.loc[index, 'stdev_cov'] = statistics.stdev(IDs), statistics.stdev(covs)
        combined_hist_scatter_plt(IDs, covs, residues_cov_x, residues_cov_y, Neff, index)
    return summary

workaround = True

if workaround:
    Nseq_path = "temp_Nseq.csv"
    tags = "a3m_files.txt"
    exp_tag = "ASFV-Georgia2007_FilterExact"
    output_dir = "/Users/jacobfenster/Documents/ASFV_Postdoc/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/clustered_at_100/stats"
    with open(tags, 'r') as f:
        MSA_list = f.readline().strip().split(' ')
    Nseq_df = pd.read_csv(Nseq_path, index_col=0)
    Nseq_dict = Nseq_df['Nseq_80'].to_dict()
    summary = MSA_summary(MSA_list, Nseq_dict)
    summary.to_csv(f"{output_dir}/{exp_tag}_summary.csv") #change dir for workaround
else:
    exp_tag, Nseq_path, MSA_list = sys.argv[1], sys.argv[2], sys.argv[3:]
    Nseq_df = pd.read_csv(Nseq_path, index_col=0)
    Nseq_dict = Nseq_df['Nseq_80'].to_dict()
    summary = MSA_summary(MSA_list, Nseq_dict)
    summary.index.name = "pMSA"
    summary.to_csv(f"{exp_tag}_summary.csv")


