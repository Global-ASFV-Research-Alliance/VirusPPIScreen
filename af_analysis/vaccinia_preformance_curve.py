# This is to plot the vaccinia MINT preformance curve 
import os, sys, glob, pdb, itertools 
import pandas as pd
import numpy as np
import analysis.AFppi as ppi
import analysis.plot as plot
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

def _full_model(control_df, tr=0.8, 
                first_col='Humphreys_maxcontact', second_col='iptm+ptm'):
    """
    This sets all rows 'second_column' equal to zero if the 'first_col' is below the 'tr'
    It also sets nan values to zero and returns a copy of the dataframe 
    """
    control_df_out = control_df.copy()
    for index in control_df_out.index:
        if control_df_out.loc[index, first_col] < tr:
            control_df_out.loc[index, second_col] = 0
        if np.isnan(control_df_out.loc[index, second_col]): 
            control_df_out.loc[index, second_col] = 0
    return control_df_out


def plot_multiple_line(Xs, Ys, data_labels, xlabel, ylabel, 
                       output_dir, basename, custom_vertical=False):
    """this plots multiple lines on the same axis
    Xs, and Ys, are lists of lists, corresponding to the same index in xlabels and ylabels
    This one is custom for Figure 1B
    """
    # Define the font properties
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 12,
           }
    legend_font = {'family': 'serif', 'weight': 'normal', 'size': 8}
    plt.figure(figsize=(2.5, 2.5))  # Set the figure size
    for X, Y, label in zip(Xs, Ys, data_labels):
        plt.plot(X, Y, label=label, linewidth=0.75) 
    ax = plt.gca()
    # Set labels with the defined font properties
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)
    # Set axis limits
    plt.xlim(0,1.1)
    plt.ylim(0,0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Set tick label sizes
    plt.xticks(fontsize=font['size']-2)
    plt.yticks(fontsize=font['size']-2)
    ax.legend(prop=legend_font, frameon=False, bbox_to_anchor=(.38, 1), loc='upper left')
    if custom_vertical:
        plt.axvline(x=0.442307692307692, color='r', linestyle='--', linewidth=0.75)
        plt.axvline(x=0.153846153846153, color='k', linestyle='--', linewidth=0.75)
    ax.xaxis.set_major_locator(MultipleLocator(.2))
    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}{basename}_plt.png")
    plt.close()

def plot_multiple_line_fig1(Xs, Ys, data_labels, xlabel, ylabel, 
                            linecolors,
                       output_dir, basename, custom_vertical=False):
    """this plots multiple lines on the same axis specific to figure 1 formatting
    Xs, and Ys, are lists of lists, corresponding to the same index in xlabels and ylabels
    This one is custom for Figure 1B
    """
    # Define the font properties
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 12,
           }
    legend_font = {'family': 'serif', 'weight': 'normal', 'size': 8}
    plt.figure(figsize=(3.2, 3.2))  # Set the figure size
    for X, Y, label, color in zip(Xs, Ys, data_labels, linecolors):
        plt.plot(X, Y, label=label, color=color, linewidth=0.75) 
    ax = plt.gca()
    # Set labels with the defined font properties
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)
    # Set axis limits
    plt.xlim(0,1.1)
    plt.ylim(0,0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Set tick label sizes
    plt.xticks(fontsize=font['size']-2)
    plt.yticks(fontsize=font['size']-2)
    legend = ax.legend(prop=legend_font, frameon=False, bbox_to_anchor=(.38, 1), loc='upper left')
    for line in legend.get_lines():
        line.set_linewidth(2)
    if custom_vertical:
        plt.axvline(x=0.442307692307692, color='#5D3A9B', linestyle='--', linewidth=0.75)
        plt.axvline(x=0.153846153846153, color='k', linestyle='--', linewidth=0.75)
    ax.xaxis.set_major_locator(MultipleLocator(.2))
    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}{basename}_plt.png", transparent=True)
    plt.close()
calculate_preformance = True 

vacc_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
if calculate_preformance:
    vacc_master = pd.read_csv(f"{vacc_out}vacc_master_partial10.csv", index_col=0)
    vacc_positive = vacc_master[vacc_master['MINT_PPI'] == 1]
    vacc_neg = vacc_master[vacc_master['MINT_PPI'] == 0]

### Using all non MINT as negative controls 

    threshold_arr = np.arange(0, 1.0001, .0001)

# calculate precision recall on just Humphreys_maxcontact algorithm
    precision_list, recall_list, df = ppi.precision_recall_curve(vacc_positive, vacc_neg, 
                                                                 threshold_arr, 'Humphreys_maxcontact')
    #df.to_csv(f"{vacc_out}HAF_mixed_vacc_preformance.csv")
# now calculate HAF -> AFmult
    vacc_positive_full = _full_model(vacc_positive)
    vacc_neg_full = _full_model(vacc_neg)
    precision_list_full, recall_list_full, df_full = ppi.precision_recall_curve(vacc_positive_full, vacc_neg_full, 
                                                                                threshold_arr, 'iptm+ptm')
    #df_full.to_csv(f"{vacc_out}full_mixed_vacc_preformance.csv")
# calculate precesion recall with RF2t++
    threshold_arr = np.arange(0, 2.001, .001)
    vacc_positive_het = vacc_master[(vacc_master['MINT_PPI'] == 1) & (vacc_master['homodimer'] == 0)]
    vacc_neg_het = vacc_master[(vacc_master['MINT_PPI'] == 0) & (vacc_master['homodimer'] == 0)]
    precision_list_rf, recall_list_rf, df_rf = ppi.precision_recall_curve(vacc_positive_het, vacc_neg_het, 
                                                                                threshold_arr, 'RF2++_paired-pMSA')
    #df_rf.to_csv(f"{vacc_out}RF2++_vacc_preformance.csv")

# try to plot by accending value 
    df_sort = df.sort_values(by='recall')
    df_full_sort = df_full.sort_values(by='recall')
    df_rf_sort = df_rf.sort_values(by='recall')
    recall_list_rf, precision_list_rf = df_rf_sort['recall'].tolist(), df_rf_sort['precision'].tolist() 
    recall_list, precision_list = df_sort['recall'].tolist(), df_sort['precision'].tolist() 
    recall_list_full, precision_list_full = df_full_sort['recall'].tolist(), df_full_sort['precision'].tolist() 

# plot preformance curve for RF2t++, mixed HAF vs full model
    Xs = [recall_list_rf, recall_list, recall_list_full]
    Ys = [precision_list_rf, precision_list, precision_list_full]
    data_labels = ['RF2t++', 'HAF', 'HAF->AF']
    xlabel, ylabel = 'recall', 'precision'
    linecolors = ['#000000', '#FFC107', '#D81B60']
    plot_multiple_line_fig1(Xs, Ys, data_labels, xlabel, ylabel, linecolors, vacc_out, 'vacc_mixed_preformance_sort_v2', custom_vertical=True)

tr_sweep = False
if tr_sweep:
### adding the line plots where the X axis is threshold!
# looking at full model only
    df = pd.read_csv(f"{vacc_out}HAF_mixed_vacc_preformance.csv", index_col=0)
    df_full = pd.read_csv(f"{vacc_out}full_mixed_vacc_preformance.csv", index_col=0)
    Xs = [df.index.tolist(), df.index.tolist()]
    Ys = [df['precision'].tolist(), df['recall'].tolist()]
    data_labels = ['precision', 'recall']
    xlabel, ylabel = 'HAF threshold', 'precision or recall'
    plot.plot_multiple_line(Xs, Ys, data_labels, xlabel, ylabel, 
                            output_dir=vacc_out, basename='HAF_trsweep_vacc_mixed',
                            xlim=(0,1), ylim=(0,1), figsize=(3, 3), linewidth=1, 
                            legendloc=(0.6, 1), x_tick=0.2, colors=['black', 'red'])

    Xs = [df_full.index.tolist(), df_full.index.tolist()]
    Ys = [df_full['precision'].tolist(), df_full['recall'].tolist()]
    data_labels = ['precision', 'recall']
    xlabel, ylabel = 'HAF->AF threshold', 'precision or recall'
    plot.plot_multiple_line(Xs, Ys, data_labels, xlabel, ylabel, 
                            output_dir=vacc_out, basename='full_trsweep_vacc_mixed',
                            xlim=(0,1), ylim=(0,1), figsize=(3, 3), linewidth=1, 
                            legendloc=(0.6, 1), x_tick=0.2, colors=['black', 'red'])




## Below are previous analysis that I am likely not going to use
if False:
    vacc_positive_full_hightr = _full_model(vacc_positive, tr=0.9)
    vacc_neg_full_hightr = _full_model(vacc_neg, tr=0.9)
    precision_list_full_h, recall_list_full_h, df_full_h = ppi.precision_recall_curve(vacc_positive_full_hightr, vacc_neg_full_hightr, 
                                                                                threshold_arr, 'iptm+ptm')
    df_full_h.to_csv(f"{vacc_out}full_mixed_tr0.9_preformance.csv")
# conduct same analysis on heterodimer only set 
    vacc_positive_het = vacc_master[(vacc_master['MINT_PPI'] == 1) & (vacc_master['homodimer'] == 0)]
    vacc_neg_het = vacc_master[(vacc_master['MINT_PPI'] == 0) & (vacc_master['homodimer'] == 0)]
    precision_list_het, recall_list_het, df_het = ppi.precision_recall_curve(vacc_positive_het, vacc_neg_het, 
                                                                             threshold_arr, 'Humphreys_maxcontact')
    df_het.to_csv(f"{vacc_out}HAF_heterodimer_vacc_preformance.csv")
    vacc_pos_het_full = _full_model(vacc_positive_het)
    vacc_neg_het_full = _full_model(vacc_neg_het)
    precision_list_het_full, recall_list_het_full, df_het_full = ppi.precision_recall_curve(vacc_pos_het_full, vacc_neg_het_full, 
                                                                                             threshold_arr, 'iptm+ptm')

    df_het_full.to_csv(f"{vacc_out}full_heterodimer_vacc_preformance.csv")
    vacc_positive_ho = vacc_master[(vacc_master['MINT_PPI'] == 1) & (vacc_master['homodimer'] == 1)]
    vacc_neg_ho = vacc_master[(vacc_master['MINT_PPI'] == 0) & (vacc_master['homodimer'] == 1)]
    precision_list_ho, recall_list_ho, df_ho = ppi.precision_recall_curve(vacc_positive_ho, vacc_neg_ho, 
                                                                             threshold_arr, 'Humphreys_maxcontact')
    df_ho.to_csv(f"{vacc_out}HAF_homodimer_vacc_preformance.csv")
    vacc_pos_ho_full = _full_model(vacc_positive_ho)
    vacc_neg_ho_full = _full_model(vacc_neg_ho)
    precision_list_ho_full, recall_list_ho_full, df_ho_full = ppi.precision_recall_curve(vacc_pos_ho_full, vacc_neg_ho_full, 
                                                                                             threshold_arr, 'iptm+ptm')
    df_ho_full.to_csv(f"{vacc_out}full_homodimer_vacc_preformance.csv")

### Using all combinations of MINT gold standards and bootstrapping them as neg controls
    mint_ppis = set()
    mint_proteins = set()
    for index in vacc_positive.index:
        protA, protB = index.split('__')[0], index.split('__')[1] 
        mint_proteins.add(protA), mint_proteins.add(protB)
        mint_ppis.add(index)

    random_pairs = list(itertools.combinations(mint_proteins, 2)) 
    mint_neg = []
    for pair in random_pairs:
        protA, protB = min([pair[0], pair[1]]), max([pair[0], pair[1]])
        index = f"{protA}__{protB}"
        if index not in mint_ppis:
            mint_neg.append(index)
    for index in mint_neg:
        vacc_master.loc[index, 'mint_neg'] = 1
    mint_neg_df = vacc_master[vacc_master['mint_neg']==1]
    bootstrap_mint_neg = mint_neg_df.sample(n=23654, replace=True)
    breakpoint()
    prec_list_mint_neg, recall_list_mint_neg, df_mint_neg = ppi.precision_recall_curve(vacc_positive, bootstrap_mint_neg,
                                                                                       threshold_arr, 'Humphreys_maxcontact')
    df_mint_neg.to_csv(f'{vacc_out}HAF_mintneg_mixed_preformance.csv')
    bootstrap_mint_neg = bootstrap_mint_neg.reset_index()
    bootstrap_mint_neg_full = _full_model(bootstrap_mint_neg)
    prec_list_mint_neg_full, recall_list_mint_neg_full, df_mint_neg_full = ppi.precision_recall_curve(vacc_positive_full, 
                                                                                                      bootstrap_mint_neg_full, 
                                                                                                      threshold_arr, 
                                                                                                      'Humphreys_maxcontact')
### doing just McCraith positive radom pairing + bootstrap as negative 
    df_mint_neg_full.to_csv(f"{vacc_out}full_mintneg_mixed_preformance.csv")
    mccriath_pos = vacc_master[vacc_master['gold_std']==1]
    mc_ppis = set()
    mc_proteins = set()

    for index in mccriath_pos.index:
        protA, protB = index.split('__')[0], index.split('__')[1] 
        mc_proteins.add(protA), mc_proteins.add(protB)
        mc_ppis.add(index)

    random_pairs = list(itertools.combinations(mc_proteins, 2)) 
    mc_neg = []
    for pair in random_pairs:
        protA, protB = min([pair[0], pair[1]]), max([pair[0], pair[1]])
        index = f"{protA}__{protB}"
        if index not in mc_ppis:
            mc_neg.append(index)
    for index in mc_neg:
        vacc_master.loc[index, 'McCraith_neg'] = 1
    mc_neg_df = vacc_master[vacc_master['McCraith_neg']==1]
    bootstrap_mc_neg = mc_neg_df.sample(n=23654, replace=True)
    prec_list_mc_neg, recall_list_mc_neg, df_mc_neg = ppi.precision_recall_curve(vacc_positive, bootstrap_mc_neg,
                                                                                       threshold_arr, 'Humphreys_maxcontact')
    df_mc_neg.to_csv(f"{vacc_out}HAF_McCraithneg_mixed_preformance.csv")
    bootstrap_mc_neg = bootstrap_mc_neg.reset_index()
    bootstrap_mc_neg_full = _full_model(bootstrap_mc_neg)
    prec_list_mc_neg_full, recall_list_mc_neg_full, df_mc_neg_full = ppi.precision_recall_curve(vacc_positive_full, 
                                                                                                      bootstrap_mc_neg_full, 
                                                                                                      threshold_arr, 
                                                                                                      'Humphreys_maxcontact')
    df_mc_neg_full.to_csv(f"{vacc_out}full_McCraithneg_mixed_preformance.csv")

    breakpoint()


