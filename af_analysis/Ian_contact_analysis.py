import os, sys, glob, pdb
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def max_contact_terimini_correction(arr):
# this definition inputs a RF2t numpy array of contact probabilities and excludes the last ten amino acids
# of the first protein and the first ten amino acids of the second protein when calculating the max contact PPI
    arr[-10:, :] = 0 # set the last 10 residues on C terminus of protein 1 to zero
    arr[:, :10] = 0 # set the first 10 residues on N terminus of protein 2 to zero
    max_contact = np.max(arr)
    return max_contact

def dir_max_contact(input_dir, protein_len_df):
    npz_files = glob.glob(f"{input_dir}*.npz")
    output_df = pd.DataFrame(columns=['max_contact_prob'])
    output_df.index.name = 'PPI'
    incorrect_shape = 0
    for npz in npz_files:
        basename = os.path.basename(npz).split('.')[0]
        index = basename.split('_info')[0]
        protein1, protein2 = basename.split('__')[0], basename.split('__')[1].split('_')[0]
        len1, len2 = int(protein_len_df.loc[protein1, 'length']), int(protein_len_df.loc[protein2, 'length'])
        data = np.load(npz)
        if data['contact_prob'].size == 0:
            breakpoint()
            continue
        contact = data['contact_prob'][:len1, len1:]
        if contact.shape[0] != len1 or contact.shape[1] != len2: #seeing some that have strange values. skipping for now 
            incorrect_shape += 1
            print(f"incorrrect shape for {npz}")
            continue
        max_contact = max_contact_terimini_correction(contact)
        output_df.loc[index, 'max_contact_prob'] = max_contact
    return output_df

def calculate_precision_recall(pos_ctrl_df, neg_ctrl_df, threshold):
    """this does not control the ratio of negative to positive controls in each df"""
    pos_above_tr = len(pos_ctrl_df[pos_ctrl_df['max_contact_prob'] > threshold])
    neg_above_tr = len(neg_ctrl_df[neg_ctrl_df['max_contact_prob'] > threshold]) 
    total_above_tr = pos_above_tr + neg_above_tr
    precision = pos_above_tr / total_above_tr
    recall = pos_above_tr / len(pos_ctrl_df.index)
    return precision, recall

def precision_recall_curve(pos_df, neg_df, threshold_arr):
    precision_list, recall_list = [], []
    for tr in threshold_arr:
        precision, recall = calculate_precision_recall(pos_df, neg_df, tr)
        precision_list.append(precision)
        recall_list.append(recall)
    return precision_list, recall_list

def plot_multiple_line(Xs, Ys, data_labels, xlabel, ylabel, output_dir,  basename):
    """this plots multiple lines on the same axis
    Xs, and Ys, are lists of lists, corresponding to the same index in xlabels and ylabels
    """
    # Define the font properties
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
           }
    plt.figure(figsize=(5, 5))  # Set the figure size
    for X, Y, label in zip(Xs, Ys, data_labels):
        plt.plot(X, Y, label=label) 
    # Set labels with the defined font properties
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)
    # Set axis limits
    plt.xlim(0,1)
    plt.ylim(0,1)
    # Set tick label sizes
    plt.xticks(fontsize=font['size'])
    plt.yticks(fontsize=font['size'])
    plt.legend()
    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}{basename}_plt.png")
    plt.close()

def violin_plot(df, x='x', y='y', hue='hue', output_dir='/home/jacob.fenster/', tag='exp_tag'):
    """ use sns here instead of matplolib https://stackoverflow.com/questions/62278350/matplotlib-seaborn-violin-plots-for-different-data-sizes""" 
    sns.set(style="whitegrid")
    plt.figure(figsize=(10,6))
    ax = sns.violinplot(data=df, x=x, y=y, hue=hue)
    plt.savefig(f"{output_dir}{tag}_plt.png")

# Gold50db has the protein-length information because this isn't in the output data file names
Gold50db = "/home/jacob.fenster/scripts/alphafold/data/Ian_MSA_proteinlen.csv"
# biological ratio of negtaive to positive PPIs
ratio = 1000 
# input 
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20240108_Ian_Gold50_paired_unpaired_analysis/"
os.makedirs(output_dir, exist_ok=True)
exp_tag = "Gold50_paired_vs_unpaired_analysis"
gold_pos_paired = "/lustrefs/fadru/projects/asfv-ppi/Output/20231221_Ian2021_Gold50_paired/"
gold_pos_paired_unpaired = "/lustrefs/fadru/projects/asfv-ppi/Output/20240102_Ian2021_Gold50_paired_unpaired/"
gold_neg_paired = "/lustrefs/fadru/projects/asfv-ppi/Output/20240102_Ian2021_neg_controls/paired/" 
gold_neg_paired_unpaired = "/lustrefs/fadru/projects/asfv-ppi/Output/20240102_Ian2021_neg_controls/paired_unpaired/"
# generate dataframes of max contact probabilities
Gold50df = pd.read_csv(Gold50db, index_col=0)
pos_paired_df = dir_max_contact(gold_pos_paired, Gold50df)
neg_paired_df = dir_max_contact(gold_neg_paired, Gold50df) 
pos_paired_unpaired_df = dir_max_contact(gold_pos_paired_unpaired, Gold50df) 
neg_paired_unpaired_df = dir_max_contact(gold_neg_paired_unpaired, Gold50df) 
# bootstrap the negative control datasets because I didn't run enough
bootstrap_neg_paired = neg_paired_df.sample(n=(len(pos_paired_df.index))*ratio, replace=True)
bootstrap_neg_paired_unpaired = neg_paired_unpaired_df.sample(n=(len(pos_paired_unpaired_df.index))*ratio, replace=True)

# array for threshold values and steps to sample
threshold_arr = np.arange(0, 1, .01)
precision_paired, recall_paired = precision_recall_curve(pos_paired_df, bootstrap_neg_paired, threshold_arr)
precision_paired_unpaired, recall_paired_unpaired = precision_recall_curve(pos_paired_unpaired_df, bootstrap_neg_paired_unpaired, threshold_arr)
Xs, Ys, data_labels = [recall_paired, recall_paired_unpaired], [precision_paired, precision_paired_unpaired], ['paired', 'paired_unpaired']
plot_multiple_line(Xs, Ys, data_labels, 'recall', 'precision', output_dir, 'Ian2021-paired_unpaired_precision-recall')

# save summary df to .csv
pos_paired_df.to_csv(f"{output_dir}Gold50_pos_ctrl_paired.csv")
neg_paired_df.to_csv(f"{output_dir}Gold50_neg_ctrl_paired.csv")
pos_paired_unpaired_df.to_csv(f"{output_dir}Gold50_pos_ctrl_paired_unpaired.csv")
neg_paired_unpaired_df.to_csv(f"{output_dir}Gold50_neg_ctrl_paired_unpaired.csv")

# generate dataframe for violin plot
columns = ['Gold50 paired', 'Gold50 paired/unpaired', 'Neg ctrl paired/unpaired']
pos_paired_df['experiment'], neg_paired_df['experiment'], pos_paired_unpaired_df['experiment'], neg_paired_unpaired_df['experiment'] = \
    'paired', 'paired', 'paired+unpaired', 'paired+unpaired'
pos_paired_df['control'], neg_paired_df['control'], pos_paired_unpaired_df['control'], neg_paired_unpaired_df['control'] = \
    'positive', 'negative', 'positive', 'negative'

pos_paired_df, neg_paired_df, pos_paired_unpaired_df, neg_paired_unpaired_df = pos_paired_df.reset_index(), neg_paired_df.reset_index(), pos_paired_unpaired_df.reset_index(), neg_paired_unpaired_df.reset_index()
violin_data = pd.concat([pos_paired_df, neg_paired_df, pos_paired_unpaired_df, neg_paired_unpaired_df])
violin_plot(violin_data, x='experiment', y='max_contact_prob', hue='control', output_dir=output_dir, tag='Gold50_paired_vs_pairedunpaired')

breakpoint()
