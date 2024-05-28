import os, sys, glob, pdb, json
import pickle
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

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

def load_pickle(file_path):
    with open(file_path, 'rb') as file:
        data = pickle.load(file)
    return data

def AF_pickle_extract(input_dir, pickle_basename):
    """
    this extracts data from an AF multimer .pkl file
    input is AF experiment directory, the file basename with extention, 
    """
    data = load_pickle(f"{input_dir}{pickle_basename}")
    ptm, iptm = float(data['ptm']), float(data['iptm'])
    return ptm, iptm

def AFmult_data_extract(exp_dirs):
    """This extracts data from the experiment sub directories of the parent input_dir"""
    df = pd.DataFrame(columns=['iptm+ptm', 'iptm', 'ptm'])
    for exp in exp_dirs:
        # make sure the folder contains results
        if os.path.exists(f"{exp}ranking_debug.json"):
            with open(f"{exp}msas/chain_id_map.json", 'r') as f:
                data = json.load(f)
            protA, protB = data['A']['description'], data['B']['description']
            # extract the iptm+ptm data
            with open(f"{exp}ranking_debug.json", 'r') as f:
                data = json.load(f)
            max_model = data['order'][0]
            max_TM_score = data['iptm+ptm'][max_model]
            # extract individual pTM and ipTM from .pkl file
            data = load_pickle(f"{exp}result_{max_model}.pkl")
            ptm, iptm = float(data['ptm']), float(data['iptm'])
            df.loc[f"{protA}__{protB}", 'iptm+ptm'] = max_TM_score
            df.loc[f"{protA}__{protB}", 'ptm'] = ptm
            df.loc[f"{protA}__{protB}", 'iptm'] = iptm
        else:
            continue
    return df

def calculate_precision_recall(pos_ctrl_df, neg_ctrl_df, threshold, column):
    """this does not control the ratio of negative to positive controls in each df"""
    pos_above_tr = len(pos_ctrl_df[pos_ctrl_df[f"{column}"] > threshold])
    neg_above_tr = len(neg_ctrl_df[neg_ctrl_df[f"{column}"] > threshold]) 
    total_above_tr = pos_above_tr + neg_above_tr
    if total_above_tr == 0:
        precision, recall = 0, 0
    else:
        precision = pos_above_tr / total_above_tr
        recall = pos_above_tr / len(pos_ctrl_df.index)
    return precision, recall

def precision_recall_curve(pos_df, neg_df, threshold_arr, column):
    """
    input the positive control and negative control df, the threshold array, and the column name of the df to analyze
    """
    precision_list, recall_list = [], []
    df = pd.DataFrame(columns=['precision', 'recall'])
    df.index.name = 'threshold'
    for tr in threshold_arr:
        precision, recall = calculate_precision_recall(pos_df, neg_df, tr, column)
        precision_list.append(precision)
        recall_list.append(recall)
        df.loc[tr, 'precision'], df.loc[tr, 'recall'] = precision, recall
    return precision_list, recall_list, df

def violin_df_formatter(data, data_columns, experiments, controls, master_column='score', master_experiment='experiment', master_control='control'):
    """this takes a dataframe with multiple columns of pTM data and outputs it in a seaborn violin plot format """
    output_df = pd.DataFrame()
    for column, experiment, control in zip(data_columns, experiments, controls):
        df = pd.DataFrame()
        df[master_column] = data[column]
        df[master_experiment] = experiment
        df[master_control] = control
        df = df.reset_index()
        output_df = pd.concat([output_df, df])
    return output_df


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
    plt.ylim(0,1.1)
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
    ax = sns.violinplot(data=df, x=x, y=y, hue=hue, cut=0)
    ax = sns.stripplot(data=df, x=x, y=y, hue=hue, jitter=True, dodge=True, color='black', size=2, edgecolor='white', linewidth=0.2)
    plt.savefig(f"{output_dir}{tag}_plt.png")

def stripplot(df, x='x', y='y', hue=None, output_dir='/home/jacob.fenster/', tag='exp_tag'):
    sns.set(stype="whitegrid")
    plt.figure(figsize=(10,6))
    if hue is None:
        ax = sns.stripplot(data=df, x=x, y=y, jitter=True, dodge=True, color='black', size=2, edgecolor='white', linewidth=0.2)
    else:
        ax = sns.stripplot(data=df, x=x, y=y, hue=hue, jitter=True, dodge=True, color='black', size=2, edgecolor='white', linewidth=0.2)
    plt.savefig(f"{output_dir}{tag}_plt.png")


def custom_merge(df1, df2, index_col, data_col1, data_col2, exp_tag1, exp_tag2):
    # this def does a custom merge allowing for either order of PPI
    for index in df1.index:
        protA, protB = df1.loc[index, index_col].split('__')[0], df1.loc[index, index_col].split('__')[1]
        index2 = f"{protB}__{protA}"
        df1.loc[index, index_col+"2"] = index2
    merge = pd.DataFrame(columns=[index_col, f"{data_col1}-{exp_tag1}", f"{data_col2}-{exp_tag2}"])
    for index1 in df1.index:
        ppi1, ppi2 = df1.loc[index1, index_col], df1.loc[index1, index_col+"2"]
        for index2 in df2.index:
            if ppi1 == df2.loc[index2, index_col] or ppi2 == df2.loc[index2, index_col]:
                merge.loc[index1] = {index_col: ppi1, f"{data_col1}-{exp_tag1}": df1.loc[index1, data_col1], f"{data_col2}-{exp_tag2}": df2.loc[index2, data_col2]}
    return merge   


def custom_merge_any(df1, df2, index_col, data_col1, data_col2):
    # this def does a custom merge allowing for either order of PPI
    # df1 and df2 must have numerical indices, the same column name for PPI indicies
    # and must have unique data columns for combination

    # extract data columns
    merge = pd.DataFrame(columns=[index_col]+data_col1+data_col2)
    for index1 in df1.index:
        ppi1, ppi2 = df1.loc[index1, index_col], f"{df1.loc[index1, index_col].split('__')[1]}__{df1.loc[index1, index_col].split('__')[0]}"
        for index2 in df2.index:
            if ppi1 == df2.loc[index2, index_col] or ppi2 == df2.loc[index2, index_col]:
                data = {}
                data[index_col] = ppi1
                for column in data_col1:
                    data[column] = df1.loc[index1, column]
                for column in data_col2:
                    data[column] = df2.loc[index2, column]
                merge.loc[index1] = data
                    
    return merge

#### Ian2021 Gold50 Analysis ####

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
precision_paired, recall_paired, ian_paired_pr_df = precision_recall_curve(pos_paired_df, bootstrap_neg_paired, threshold_arr, 'max_contact_prob')
precision_paired_unpaired, recall_paired_unpaired, ian_paired_pr_df = precision_recall_curve(pos_paired_unpaired_df, bootstrap_neg_paired_unpaired, threshold_arr, 'max_contact_prob')
Xs, Ys, data_labels = [recall_paired, recall_paired_unpaired], [precision_paired, precision_paired_unpaired], ['paired', 'paired_unpaired']
plot_multiple_line(Xs, Ys, data_labels, 'recall', 'precision', output_dir, 'Ian2021-paired_unpaired_precision-recall')

# save summary df to .csv
pos_paired_df.to_csv(f"{output_dir}Gold50_pos_ctrl_paired.csv")
neg_paired_df.to_csv(f"{output_dir}Gold50_neg_ctrl_paired.csv")
pos_paired_unpaired_df.to_csv(f"{output_dir}Gold50_pos_ctrl_paired_unpaired.csv")
neg_paired_unpaired_df.to_csv(f"{output_dir}Gold50_neg_ctrl_paired_unpaired.csv")

# generate dataframe for violin plot
pos_paired_df['experiment'], neg_paired_df['experiment'], pos_paired_unpaired_df['experiment'], neg_paired_unpaired_df['experiment'] = \
    'paired', 'paired', 'paired+unpaired', 'paired+unpaired'
pos_paired_df['control'], neg_paired_df['control'], pos_paired_unpaired_df['control'], neg_paired_unpaired_df['control'] = \
    'positive', 'negative', 'positive', 'negative'

pos_paired_df, neg_paired_df, pos_paired_unpaired_df, neg_paired_unpaired_df = pos_paired_df.reset_index(), neg_paired_df.reset_index(), pos_paired_unpaired_df.reset_index(), neg_paired_unpaired_df.reset_index()
violin_data = pd.concat([pos_paired_df, neg_paired_df, pos_paired_unpaired_df, neg_paired_unpaired_df])
violin_plot(violin_data, x='experiment', y='max_contact_prob', hue='control', output_dir=output_dir, tag='Ian_Gold50_violin_paired_vs_pairedunpaired')

#### AF Gold50  multimer analysis ####

mult_Gold50_input_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-pos_ctrl/"
mult_neg_ctrl_input_dir_1 = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-neg_ctrl_short_reuse_MSA/"
mult_neg_ctrol_input_dir_2 = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-neg_ctrl_short/"

mult_Gold50_exps = glob.glob(f"{mult_Gold50_input_dir}*/")
mult_neg_exp1 = glob.glob(f"{mult_neg_ctrl_input_dir_1}*/")
mult_neg_exp2 = glob.glob(f"{mult_neg_ctrol_input_dir_2}*/")
mult_neg_exps = mult_neg_exp1 + mult_neg_exp2
#extract data
mult_Gold50_df = AFmult_data_extract(mult_Gold50_exps)    
mult_neg_df = AFmult_data_extract(mult_neg_exps)    
bootstrap_mult_neg = mult_neg_df.sample(n=(len(mult_Gold50_df.index))*ratio, replace=True)

threshold_arr = np.arange(0, 1, .01)

prec_mult_iptm, recall_mult_iptm, mult_pc_df_iptm = precision_recall_curve(mult_Gold50_df, bootstrap_mult_neg, threshold_arr, 'iptm')
prec_mult_ptm, recall_mult_ptm, mult_pc_df_ptm = precision_recall_curve(mult_Gold50_df, bootstrap_mult_neg, threshold_arr, 'ptm')
prec_mult_iptm_ptm, recall_mult_iptm_ptm, mult_pc_df_iptm_ptm = precision_recall_curve(mult_Gold50_df, bootstrap_mult_neg, threshold_arr, 'iptm+ptm')
# plot precision recall for the three TM metrics
Xs, Ys, data_labels = [recall_mult_ptm, recall_mult_iptm, recall_mult_iptm_ptm], [prec_mult_ptm, prec_mult_iptm, prec_mult_iptm_ptm], ['pTM', 'ipTM', 'ipTM+pTM']
plot_multiple_line(Xs, Ys, data_labels, 'recall', 'precision', output_dir, 'AFmult_TMcompare_precisionrecall')
#save data
mult_Gold50_df.to_csv(f"{output_dir}AFmult_Gold50_TMdata.csv")
mult_neg_df.to_csv(f"{output_dir}AFmult_neg_TMdata.csv")
#plot violin of AF mult data

Gold50_violin_data = violin_df_formatter(mult_Gold50_df, ['iptm', 'ptm', 'iptm+ptm'], ['iptm', 'ptm', 'iptm+ptm'], ['positive', 'positive', 'positive'])
neg_violin_data = violin_df_formatter(mult_neg_df, ['iptm', 'ptm', 'iptm+ptm'], ['iptm', 'ptm', 'iptm+ptm'], ['negative', 'negative', 'negative'])
AF_violin_data = pd.concat([Gold50_violin_data, neg_violin_data])
violin_plot(AF_violin_data, x='experiment', y='score', hue='control', output_dir=output_dir, tag='AFmult_TMpanelscores_violin')
AF_violin_data.to_csv(f"{output_dir}AFmult_violindata.csv")
#### look at threshold vs precision/recall plots to pick tr ####
# plot Ian paired
Xs, Ys, data_labels = [threshold_arr, threshold_arr], [recall_paired, precision_paired], ['recall Ip', 'precision Ip']
plot_multiple_line(Xs, Ys, data_labels, 'threshold', 'score', output_dir, 'Ian_Paired_pc_vs_threshold')
# plot Ian paired+unpaired
Xs, Ys, data_labels = [threshold_arr, threshold_arr], [recall_paired_unpaired, precision_paired_unpaired], ['recall Ipu', 'precision Ipu']
plot_multiple_line(Xs, Ys, data_labels, 'threshold', 'score', output_dir, 'Ian_Paired-unpaired_pc_vs_threshold')
# plot AFmult iptm+ptm
Xs, Ys, data_labels = [threshold_arr, threshold_arr], [recall_mult_iptm_ptm, prec_mult_iptm_ptm], ['recall AFmult i+p', 'precision AFmult i+p']
plot_multiple_line(Xs, Ys, data_labels, 'threshold', 'score', output_dir, 'AFmult_pc_vs_threshold')

#### Union analysis to compare precision recall curves across experiments ####
pos_paired_df.rename(columns={'max_contact_prob':'max_contact_prob-paired'}, inplace=True), neg_paired_df.rename(columns={'max_contact_prob':'max_contact_prob-paired'}, inplace=True)
pos_paired_unpaired_df.rename(columns={'max_contact_prob':'max_contact_prob-paired_unpaired'}, inplace=True), neg_paired_unpaired_df.rename(columns={'max_contact_prob':'max_contact_prob-paired_unpaired'}, inplace=True)
Ian_pos_merge = custom_merge_any(pos_paired_df, pos_paired_unpaired_df, 'PPI', ['max_contact_prob-paired', 'control'], ['max_contact_prob-paired_unpaired'])
Ian_neg_merge = custom_merge_any(neg_paired_df, neg_paired_unpaired_df, 'PPI', ['max_contact_prob-paired', 'control'], ['max_contact_prob-paired_unpaired'])
# bootstrap the neg union
Ian_neg_merge_boot = Ian_neg_merge.sample(n=(len(Ian_pos_merge.index))*ratio, replace=True)
#plot Ian2021 union Precision recall plot
prec_ian_union_p, rec_ian_union_p, ian_union_pr_p_df = precision_recall_curve(Ian_pos_merge, Ian_neg_merge_boot, threshold_arr, 'max_contact_prob-paired')
prec_ian_union_pu, rec_ian_union_pu, ian_union_pr_pu_df = precision_recall_curve(Ian_pos_merge, Ian_neg_merge_boot, threshold_arr, 'max_contact_prob-paired_unpaired')
#plot the comparason of union
Xs, Ys, data_labels = [rec_ian_union_p, rec_ian_union_pu], [prec_ian_union_p, prec_ian_union_pu], ['paired', 'paired_unpaired']
plot_multiple_line(Xs, Ys, data_labels, 'recall', 'precision', output_dir, 'Ian2021-Union-paired_unpaired_precision-recall')

mult_Gold50_df, mult_neg_df = mult_Gold50_df.reset_index(), mult_neg_df.reset_index()
mult_Gold50_df.rename(columns={'index':'PPI'}, inplace=True), mult_neg_df.rename(columns={'index':'PPI'}, inplace=True)
pos_union = custom_merge_any(Ian_pos_merge, mult_Gold50_df, 'PPI', ['max_contact_prob-paired', 'max_contact_prob-paired_unpaired'], ['iptm+ptm', 'iptm', 'ptm'])
neg_union = custom_merge_any(Ian_neg_merge, mult_neg_df, 'PPI', ['max_contact_prob-paired', 'max_contact_prob-paired_unpaired'], ['iptm+ptm', 'iptm', 'ptm'])
neg_union_boot = neg_union.sample(n=(len(pos_union.index))*ratio, replace=True)
#All union PR curves 
prec_all_union_ian_p, rec_all_union_ian_p, all_union_ian_p_df = precision_recall_curve(pos_union, neg_union_boot, threshold_arr, 'max_contact_prob-paired')
prec_all_union_ian_pu, rec_all_union_ian_pu, all_union_ian_pu_df = precision_recall_curve(pos_union, neg_union_boot, threshold_arr, 'max_contact_prob-paired_unpaired')
prec_all_union_mult, rec_all_union_mult, all_union_mult_df = precision_recall_curve(pos_union, neg_union_boot, threshold_arr, 'iptm+ptm')

Xs, Ys, data_labels = [rec_all_union_ian_p, rec_all_union_ian_pu, rec_all_union_mult], [prec_all_union_ian_p, prec_all_union_ian_pu, prec_all_union_mult], ['paired', 'paired_unpaired', 'AFmult iptm+ptm']
plot_multiple_line(Xs, Ys, data_labels, 'recall', 'precision', output_dir, 'Union_all-recall')

# looking at high score Gold50 neg
AFmult_high_neg_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/high_score_neg/"
AFmult_high_neg_dirs = glob.glob(f"{AFmult_high_neg_dir}*/")
mult_high_neg = AFmult_data_extract(AFmult_high_neg_dirs)
mult_high_neg.index.name = "PPI"
mult_high_neg = mult_high_neg.reset_index()
hineg_union = custom_merge_any(mult_high_neg, neg_paired_unpaired_df, 'PPI', ['iptm+ptm', 'iptm', 'ptm'], ['max_contact_prob-paired_unpaired'])
hineg_union = custom_merge_any(hineg_union, neg_paired_df, 'PPI', ['iptm+ptm', 'iptm', 'ptm', 'max_contact_prob-paired_unpaired'], ['max_contact_prob-paired'])
hineg_union.to_csv(f"{output_dir}high_score_neg_union.csv")
breakpoint()
hineg_violin = violin_df_formatter(hineg_union, ['max_contact_prob-paired', 'max_contact_prob-paired_unpaired', 'iptm+ptm'], ['max_contact_prob-paired', 'max_contact_prob-paired_unpaired', 'iptm+ptm'], ['negative', 'negative', 'negative'])
stripplot(hineg_union, x='score', y='experiment', hue=None, output_dir=output_dir, tag='high_score_neg_stripplot')

