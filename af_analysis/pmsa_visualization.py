import os, sys, pdb
import pandas as pd
import numpy as np
import analysis.plot as plot 

vacc_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
asfv_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
vacc_master = pd.read_csv(f"{vacc_out}vacc_master_partial10.csv")
asfv_master = pd.read_csv(f"{asfv_out}asfv_master_partial9.csv")

### Heterodimer work

# Vaccinia work
# plot histograms of Neff80 paired 
vacc_master_het = vacc_master[vacc_master['homodimer']==0]
plot.plot_single_hist(vacc_master_het, 'Neff80-paired', f'{vacc_out}vacc_het_Neff80paired_hist.png', 
                         xlim=(0,20), xlabel='Neff80-paired', n_bins=800, title=None)
# plot Nseq80 paired 
vacc_master_het = vacc_master[vacc_master['homodimer']==0]
plot.plot_single_hist(vacc_master_het, 'Nseq80-paired', f'{vacc_out}vacc_het_Nseq80paired_hist.png', 
                         xlim=(0, 250), xlabel='Nseq80-paired', n_bins=1300, title=None, yaxis='log')
# plot hist of Nseq paired
plot.plot_single_hist(vacc_master_het, 'Nseq-paired', f'{vacc_out}vacc_het_Nseqpaired_hist.png', 
                         xlim=(0, 500), xlabel='Nseq-paired', n_bins=800, title=None, yaxis='log')
# ASFV work
asfv_master_het = asfv_master[asfv_master['homodimer']==0]
plot.plot_single_hist(asfv_master_het, 'Neff80-paired', f'{asfv_out}asfv_het_Neff80paired_hist.png', 
                         xlim=(0,20), xlabel='Neff80-paired', n_bins=50, title=None)
# plot Nseq80 paired 
asfv_master_het = asfv_master[asfv_master['homodimer']==0]
plot.plot_single_hist(asfv_master_het, 'Nseq80-paired', f'{asfv_out}asfv_het_Nseq80paired_hist.png', 
                         xlim=(0,250), xlabel='Nseq80-paired', n_bins=100, title=None, yaxis='log')
# plot hist of Nseq paired
plot.plot_single_hist(asfv_master_het, 'Nseq-paired', f'{asfv_out}asfv_het_Nseqpaired_hist.png', 
                         xlim=(0, 500), xlabel='Nseq-paired', n_bins=100, title=None, yaxis='log')

### Homodimer work

# Vaccinia work
# plot histograms of protA unpaired (as they are mirror images)
vacc_master_ho = vacc_master[vacc_master['homodimer']==1]
plot.plot_single_hist(vacc_master_ho, 'Neff80-protA', f'{vacc_out}vacc_ho_Neff80_unpaired_hist.png', 
                      xlim=(0, 700), xlabel='Neff80-unpaired', n_bins=20)

plot.plot_single_hist(vacc_master_ho, 'Nseq-protA', f'{vacc_out}vacc_ho_Nseq_unpaired_hist.png', 
                      xlim=(0, 11000), xlabel='Nseq-unpaired', n_bins=20)
# ASFV work
asfv_master_ho = asfv_master[asfv_master['homodimer']==1]
plot.plot_single_hist(asfv_master_ho, 'Neff80-protA', f'{asfv_out}asfv_ho_Neff80_unpaired_hist.png', 
                      xlim=(0, 700), xlabel='Neff80-unpaired', n_bins=20)

plot.plot_single_hist(asfv_master_ho, 'Nseq-protA', f'{asfv_out}asfv_ho_Nseq_unpaired_hist.png', 
                      xlim=(0, 11000), xlabel='Nseq-unpaired', n_bins=20)

