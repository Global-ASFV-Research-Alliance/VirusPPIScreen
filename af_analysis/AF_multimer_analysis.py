import os, sys, glob, pdb, json
import pickle
import numpy as np
import pandas as pd 
import seaborn as sns

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
    precision = pos_above_tr / total_above_tr
    recall = pos_above_tr / len(pos_ctrl_df.index)
    return precision, recall

def precision_recall_curve(pos_df, neg_df, threshold_arr, column):
    """
    input the positive control and negative control df, the threshold array, and the column name of the df to analyze
    """
    precision_list, recall_list = [], []
    for tr in threshold_arr:
        precision, recall = calculate_precision_recall(pos_df, neg_df, tr, column)
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
    plot.legend()
    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}{basename}_plt.png")
    plt.close()

mult_Gold50_input_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-pos_ctrl/"
mult_neg_ctrl_input_dir_1 = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-neg_ctrl_short_reuse_MSA/"
mult_neg_ctrol_input_dir_2 = "/lustrefs/fadru/projects/asfv-ppi/Output/20231222_AFmult_Gold50/Gold50-neg_ctrl_short/"

mult_Gold50_exps = glob.glob(f"{mult_Gold50_input_dir}*/")
mult_neg_exp1 = glob.glob(f"{mult_neg_ctrl_input_dir_1}*/")
mult_neg_exp2 = glob.glob(f"{mult_neg_ctrol_input_dir_2}*/")
mult_neg_exps = mult_neg_exp1 + mult_neg_exp2
mult_Gold50_df = AFmult_data_extract(mult_Gold50_exps)    
mult_neg_df = AFmult_data_extract(mult_neg_exps)    
breakpoint()
bootstrap_mult_neg = neg_df.sample(n=(len(mult_Gold50_df.index))*ratio, replace=True)

threshold_arr = np.arange(0, 1, .01)
