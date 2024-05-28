# This script is to visualize the master data tables from the ASFV, Vaccinia, and FMDV AI protein prediction work
import os, sys, glob, pdb
import numpy as np
import pandas as pd
import analysis.plot as plot 
import analysis.AFppi as ppi


### DATA INPUT ###
vacc_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
fmdv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV_all_homo-heterodimer/"
asfv_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
combined_analysis_out = "/lustrefs/fadru/projects/asfv-ppi/data/combined_virus_analysis/"
# data column labels 
ian2021_col = 'Humphreys_maxcontact' # name for Ian Humphrey's 2021 max contact calcuation
afmult_col = 'iptm+ptm' # name for AFmultimer iptm+ptm calculation
rf_mega_col = 'RF2++_mega-pMSA'
rf_paired_col = 'RF2++_paired-pMSA'

columns = ['PPI_A', 'PPI_B', 'proteinA', 'proteinB', ian2021_col, afmult_col, 
           'gold_std', 'seqA', 'seqB', 'protA_description', 'protB_description',
           'Neff80-paired', 'Neff80-protA', 'Neff80-protB']

def _label_homodimers(master):
    for index in master.index:
        protA, protB = index.split('__')[0], index.split('__')[1]
        if protA == protB:
            master.loc[index, 'homodimer'] = 1 
        else:
            master.loc[index, 'homodimer'] = 0 
    return master

def _master_heterodimer_summary_metrics(master, row_label='test', homodimer=0):
    summary_dict = {'total_ppi': (master['homodimer']==homodimer).sum(), 
                    f'{ian2021_col}>=0.8': len(master[(master[ian2021_col] >= 0.8) & (master['homodimer']==homodimer)]), 
                    f'{ian2021_col}>=0.9': len(master[(master[ian2021_col] >= 0.9) & (master['homodimer']==homodimer)]), 
                    f'AFmult_total': ((master[afmult_col].notna()) & (master['homodimer']==homodimer)).sum(),
                    f'AFmult-{afmult_col}>=0.8': len(master[(master[afmult_col] >= 0.8) & (master['homodimer']==homodimer)]), 
                    f'AFmult-{afmult_col}>=0.9': len(master[(master[afmult_col] >= 0.9) & (master['homodimer']==homodimer)])}
    df = pd.DataFrame(data=summary_dict, index=[row_label])
    return df
                    
def _alphabetize_ppi(df):
    # this parses the ppi protA and protB from the index col and renames according to alphabetical order
    for index in df.index:
        prot1, prot2 = index.split('__')[0], index.split('__')[1]
        protA, protB = min([prot1, prot2]), max([prot1, prot2])
        new_index = f"{protA}__{protB}"
        df.rename(index={index: new_index}, inplace=True)
    return df

### ASFV work ###
if False:
    asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial8.csv", index_col=0)
# compare RF to AF 
    rf_mega, rf_paired = asfv_master[rf_mega_col].tolist(), asfv_master[rf_paired_col].tolist()
    ian2021 = asfv_master[ian2021_col].tolist()
    x_label, y_label = rf_mega_col, ian2021_col
    out_file = f"{asfv_analysis_out}asfv_Ian2021_vs_rf-mega_scatter.png"
    plot.plot_scatter_linear_reg(rf_mega, ian2021, out_file, x_label, y_label, plot_type='ppt')

    x_label, y_label = rf_paired_col, ian2021_col
    out_file = f"{asfv_analysis_out}asfv_Ian2021_vs_rf-paired_scatter.png"
    plot.plot_scatter_linear_reg(rf_paired, ian2021, out_file, x_label, y_label, plot_type='ppt')
# calculate summary metrics on heterodimer set
    asfv_master = _label_homodimers(asfv_master)
    asfv_master.to_csv(f"{asfv_analysis_out}asfv_master_partial9.csv")
    asfv_het_summary = _master_heterodimer_summary_metrics(asfv_master, 'ASFV', homodimer=0)
    asfv_ho_summary = _master_heterodimer_summary_metrics(asfv_master, 'ASFV', homodimer=1)
# ASFV gold standard work
    asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial9.csv", index_col=0)
    asfv_gold_paired = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/ASFV/early_analysis/ASFV_expGold_Ian2021_paired_data.csv", index_col='PPI')
    asfv_gold_paired = asfv_gold_paired.drop('Unnamed: 0', axis=1)
    asfv_gold_paired = _alphabetize_ppi(asfv_gold_paired)
    asfv_gold_master = pd.DataFrame(columns=[f'{ian2021_col}-mega', f'{ian2021_col}-paired', afmult_col])

    asfv_gold_master = ppi._load_master_index(asfv_master[asfv_master['gold_std']==1], asfv_gold_master, [afmult_col], [afmult_col])
    asfv_gold_master = ppi._load_master_index(asfv_master[asfv_master['gold_std']==1], asfv_gold_master, [ian2021_col], [f'{ian2021_col}-mega'])
    asfv_gold_master = ppi._load_master_index(asfv_gold_paired, asfv_gold_master, ['max_contact_prob'], [f'{ian2021_col}-paired'])
    asfv_gold_master.to_csv(f"{asfv_analysis_out}asfv_gold_std_master.csv")

    conditions=[{'col_name': f'{ian2021_col}-paired', 'fig_name': 'HAF-only_paired', 'control': 'pos'},
                {'col_name': f'{ian2021_col}-mega', 'fig_name': 'HAF-mega', 'control': 'pos'},
                {'col_name': afmult_col, 'fig_name': 'AFmult', 'control': 'pos'}]
    asfv_gold_violin = plot.violin_master_formatter(asfv_gold_master, conditions)
    plot.violin_plot(asfv_gold_violin, output_file=f'{asfv_analysis_out}asfv_gold_compare_violin.png', hue=None)
# Nseq80 work
asfv_master = pd.read_csv(f"{asfv_analysis_out}asfv_master_partial9.csv", index_col=0)
asfv_gold_het = asfv_master[(asfv_master['gold_std'] == 1) & (asfv_master['homodimer'] == 0)]
asfv_gold_mega = asfv_gold_het[ian2021_col].tolist()
asfv_gold_nseq_paired = asfv_gold_het['Nseq80-paired'].tolist()
asfv_gold_nseq_protA = asfv_gold_het['Nseq80-protA'].tolist()
asfv_gold_nseq_protB = asfv_gold_het['Nseq80-protB'].tolist()
plot.plot_scatter_linear_reg(asfv_gold_nseq_paired, asfv_gold_mega, f"{asfv_analysis_out}Ian_mega-vs-Neffpaired.png",
                             x_label='log(Neff80-paired)', y_label='HAF-mega', xlim=(10**0, 10**4), 
                             tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
plot.plot_scatter_linear_reg(asfv_gold_nseq_protA, asfv_gold_mega, f"{asfv_analysis_out}Ian_mega-vs-NeffprotA.png",
                             x_label='log(Neff80-protA)', y_label='HAF-mega', xlim=(10**0, 10**4), 
                             tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
plot.plot_scatter_linear_reg(asfv_gold_nseq_protB, asfv_gold_mega, f"{asfv_analysis_out}Ian_mega-vs-NeffprotB.png",
                             x_label='log(Neff80-protB)', y_label='HAF-mega', xlim=(10**0, 10**4), 
                             tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')



### Vaccinia Work ###
if False:
    vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial7.csv", index_col=0)
# compare RF to AF
    rf_mega, rf_paired = vacc_master[rf_mega_col].tolist(), vacc_master[rf_paired_col].tolist()
    ian2021 = vacc_master[ian2021_col].tolist()
    x_label, y_label = rf_mega_col, ian2021_col
    out_file = f"{vacc_analysis_out}vacc_Ian2021_vs_rf-mega_scatter.png"
    plot.plot_scatter_linear_reg(rf_mega, ian2021, out_file, x_label, y_label, plot_type='ppt')

    x_label, y_label = rf_paired_col, ian2021_col
    out_file = f"{vacc_analysis_out}vacc_Ian2021_vs_rf-paired_scatter.png"
    plot.plot_scatter_linear_reg(rf_paired, ian2021, out_file, x_label, y_label, plot_type='ppt')
    vacc_master = _label_homodimers(vacc_master)
    vacc_master.to_csv(f"{vacc_analysis_out}vacc_master_partial8.csv")
    vacc_het_summary = _master_heterodimer_summary_metrics(vacc_master, 'Vaccinia', homodimer=0)
    vacc_ho_summary = _master_heterodimer_summary_metrics(vacc_master, 'Vaccinia', homodimer=1)
# vaccinia Gold standard work
    mccraith_afmult = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_AFmult.csv", index_col=0)
    mccraith_het = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_Ian_heterodimer_paired-unpaired.csv", index_col=0)
    mccraith_het_paired = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_Ian_heterodimer_paired.csv", index_col=0)
    mccraith_ho = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000/McCraith2000_Ian_homodimer.csv", index_col=0)
#clean up afmult homodimer names
    for index in mccraith_afmult.index:
        if len(index.split('__')) == 4:
            protA, protB = index.split('__')[0], index.split('__')[2]
            new_index = f"{protA}__{protB}"
            mccraith_afmult.rename(index={index: new_index}, inplace=True)
# alphabetize PPIs for normal indexing
    mccraith_afmult = _alphabetize_ppi(mccraith_afmult)
    mccraith_het = _alphabetize_ppi(mccraith_het)
    mccraith_ho = _alphabetize_ppi(mccraith_ho)
    mccraith_het_paired = _alphabetize_ppi(mccraith_het_paired)
    mccraith_master = pd.DataFrame(columns=[f'{ian2021_col}-mega', f'{ian2021_col}-paired', afmult_col])
    mccraith_master = ppi._load_master_index(mccraith_afmult, mccraith_master, [afmult_col], [afmult_col])
    mccraith_master = ppi._load_master_index(mccraith_het, mccraith_master, ['max_contact_prob'], [f'{ian2021_col}-mega'])
    mccraith_master = ppi._load_master_index(mccraith_het_paired, mccraith_master, ['max_contact_prob'], [f'{ian2021_col}-paired'])
    mccraith_master = ppi._load_master_index(mccraith_ho, mccraith_master, ['max_contact_prob'], [f'{ian2021_col}-mega'])
    mccraith_master = _label_homodimers(mccraith_master)
    mccraith_master.to_csv(f"{vacc_analysis_out}McCraith2000_master.csv")
# scatter of mega paired vs mega paired-unpaired 
    mccraith_master_het = mccraith_master[mccraith_master['homodimer']==0]
    mega = mccraith_master_het[f'{ian2021_col}-mega'].tolist()
    paired = mccraith_master_het[f'{ian2021_col}-paired'].tolist()
    x_label, y_label = f'{ian2021_col}-paired', f'{ian2021_col}-mega' 
    out_file = f"{vacc_analysis_out}McCraith_paired_vs_mega-pu_scatter.png"
    plot.plot_scatter_linear_reg(paired, mega, out_file, x_label, y_label, plot_type='ppt')
# seaborn plot with individual datapoints for clarity
    conditions=[{'col_name': f'{ian2021_col}-paired', 'fig_name': 'HAF-only_paired', 'control': 'pos'},
                {'col_name': f'{ian2021_col}-mega', 'fig_name': 'HAF-mega', 'control': 'pos'},
                {'col_name': afmult_col, 'fig_name': 'AFmult', 'control': 'pos'}]
    mccraith_violin = plot.violin_master_formatter(mccraith_master_het, conditions)
    plot.violin_plot(mccraith_violin, output_file=f'{vacc_analysis_out}McCraith_compare_violin.png', hue=None)

    vacc_master = pd.read_csv(f"{vacc_analysis_out}vacc_master_partial8.csv", index_col=0)
# plot neff vs mega and AFmult for gold standard set 
    mccraith_master_het = vacc_master[(vacc_master['gold_std']==1) & (vacc_master['homodimer']==0)]
    mccraith_het_afmult = mccraith_master_het[afmult_col].tolist()
    mccraith_het_mega = mccraith_master_het[ian2021_col].tolist()
    mccraith_neff_paired = mccraith_master_het['Neff80-paired'].tolist()
    mccraith_neff_protA = mccraith_master_het['Neff80-protA'].tolist()
    mccraith_neff_protB = mccraith_master_het['Neff80-protB'].tolist()
    plot.plot_scatter_linear_reg(mccraith_neff_paired, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-Neffpaired.png",
                                 x_label='log(Neff80-paired)', y_label='HAF-mega', xlim=(10**-2, 3*10**2), 
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
    plot.plot_scatter_linear_reg(mccraith_neff_protA, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-NeffprotA.png",
                                 x_label='log(Neff80-protA)', y_label='HAF-mega', xlim=(10**-2,3*10**2),
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
    plot.plot_scatter_linear_reg(mccraith_neff_protB, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-NeffprotB.png",
                                 x_label='log(Neff80-protB)', y_label='HAF-mega', xlim=(10**-2,3*10**2),
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
# Plot Nseq80
    mccraith_neff_paired = mccraith_master_het['Nseq80-paired'].tolist()
    mccraith_neff_protA = mccraith_master_het['Nseq80-protA'].tolist()
    mccraith_neff_protB = mccraith_master_het['Nseq80-protB'].tolist()
    plot.plot_scatter_linear_reg(mccraith_neff_paired, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-Nseqpaired.png",
                                 x_label='log(Nseq80-paired)', y_label='HAF-mega', xlim=(10**0, 10**4),
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
    plot.plot_scatter_linear_reg(mccraith_neff_protA, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-NseqprotA.png",
                                 x_label='log(Nseq80-protA)', y_label='HAF-mega', xlim=(10**0, 10**4), 
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')
    plot.plot_scatter_linear_reg(mccraith_neff_protB, mccraith_het_mega, f"{vacc_analysis_out}Ian_mega-vs-NseqprotB.png",
                                 x_label='log(Nseq80-protB)', y_label='HAF-mega', xlim=(10**0, 10**4), 
                                 tick_interval=None, xscale='log', lin_reg=False, plot_type='ppt_smalldata')

breakpoint()
### combining viruses 
if False:
    combined_het_summary = pd.concat([vacc_het_summary, asfv_het_summary], ignore_index=False)
    combined_ho_summary = pd.concat([vacc_ho_summary, asfv_ho_summary], ignore_index=False)
    combined_het_summary.to_csv(f"{combined_analysis_out}heterodimer_summary_metrics.csv")
    combined_ho_summary.to_csv(f"{combined_analysis_out}homodimer_summary_metrics.csv")



### pMSA analysis. Length, size, Neff. 

### Precision recall curves of the yeast work
# might need more data for this analysis

### Gold standard analysis. ASFV, McCraith, FMDV Y2H
# bar plots of percent significant calls? with Ian2021 vs AFmult


### Number of proteins and their names above thresholds
