import sys, os, pdb
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def create_annot(ax, sc, labels, ind):
    pos = sc.get_offsets()[ind["ind"][0]]
    annot = ax.annotate("", xy=pos, xytext=(10,10), textcoords="offset points",
                        arrowprops=dict(arrowstyle="->", lw=0.5), fontsize=8)
    text = ', '.join([labels[n] for n in ind['ind']])
    annot.set_text(text)
    return annot

def onclick(event, ax, sc, annot_list, labels):
    if event.inaxes == ax:
        cont, ind = sc.contains(event)
        if cont:
            clicked_pos = sc.get_offsets()[ind["ind"][0]]

            # Check if an annotation for this point already exists
            for annot in annot_list:
                # Compare each element of the positions
                if all(annot.xy == clicked_pos):
                    annot.set_visible(False)  # Hide the annotation
                    annot_list.remove(annot)  # Remove it from the list
                    plt.draw()
                    return

            # If not, create a new annotation
            annot = create_annot(ax, sc, labels, ind)
            annot.set_visible(True)
            annot_list.append(annot)
            plt.draw()

def plot_scatter(X, Y, x_label, y_label, xlim, ylim, custom_labels):
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 18,
           }

    fig, ax = plt.subplots(figsize=(5, 5))  # Set the figure size
    sc = ax.scatter(X, Y, s=2)  # Create a scatter plot
    ax.set_xlabel(x_label, fontdict=font)  # Set the x-axis label
    ax.set_ylabel(y_label, fontdict=font)  # Set the y-axis label
    ax.tick_params(axis='both', labelsize=font['size'])  # Set tick label sizes
    plt.tight_layout()  # Adjust the layout to prevent clipping
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    annot_list = []
    fig.canvas.mpl_connect("button_press_event", lambda event: onclick(event, ax, sc, annot_list, custom_labels))

    plt.show()

# Data input
results_file = "/Users/jacobfenster/Documents/ASFV_Postdoc/NPZ_analysis/output/jackhmmer_all_virus/ASFV_ASFV_clust99/ASFV_ASFV_clust99_PPIscores_combined.csv"
results_df = pd.read_csv(results_file, index_col=0)
#extract data labels
labels = []
for index in results_df.index:
    labels.append(f"{index.split('-AA')[0]}::{index.split('__')[1].split('-AA')[0]}")
#plot signal vs various parameters
RF2 = results_df['RF2t+ max contact probability'].tolist()
RF2plus = results_df['RF2t++ scores'].tolist()
Neff = results_df['Neff_80'].tolist()
Naln = results_df['#_alns'].tolist()
Nclust = results_df['Nseq_80'].tolist()
#plot_scatter(Neff, RF2plus, 'Neff_90', 'RF2+ max contact prob', 'ASFV_clust99_90', 'Neff_vs_RF2')
#plot_scatter(Naln, RF2, '# aln', 'RF2+ max contact prob', 'ASFV_clust99_90', 'NumAln_vs_RF2')
xlim_clust = (0, 50)
ylim_RF2 = (0.85, max(RF2)+.1)
plot_scatter(Nclust, RF2, 'N_clust 80', 'RF2+ max contact prob', xlim_clust, ylim_RF2, labels)
ylim_RF2plus = (0, max(RF2plus)+.1)
#plot_scatter(Nclust, RF2plus, 'N_clust 80', 'RF2+ max contact prob', xlim_clust, ylim_RF2plus, labels)

