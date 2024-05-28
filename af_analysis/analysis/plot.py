import os, sys, glob, json, pdb
import pickle
import numpy as np
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

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

def violin_master_formatter(master,
                     conditions=[{'col_name': 'name', 
                                 'fig_name': 'name', 
                                 'control': 'pos/neg'}]): 
    """
    this formats a pandas df into a format for seaborn violin plot visualization. 
    conditions is a list of dictionaries where the keys are below
    'col_name': name of the column in the input df
    'fig_name': name displayed on the violin plot x axis
    'control': if the given condition is a positive or negative control. can set all to same value and ignore later
    outputs a df with numerical indicies, a PPI column, a score column with the col_name data, 
    an experiment column with the 'fig_name' input, and a control column with 'control' data 
    master must contain all col_names and have the index col be the PPI 
    """
    violin_df = pd.DataFrame(columns=['PPI', 'score', 'experiment', 'control'])
    for condition in conditions:
        condition_data = {'PPI': master.index.tolist(),
                          'score': master[condition['col_name']].tolist(),
                          'experiment': [condition['fig_name']] * len(master.index) ,
                          'control': [condition['control']] * len(master.index)}
        df = pd.DataFrame(data=condition_data)
        violin_df = pd.concat([violin_df, df], ignore_index=True) 
    return violin_df 

def plot_multiple_line(Xs, Ys, data_labels, xlabel, ylabel, 
                       output_dir, basename,
                       xlim=(0,1), ylim=(0,1), figsize=(2.5, 2.5), linewidth=0.75, 
                       legendloc=(0.4, 1), x_tick=0.2, colors=None):
    """this plots multiple lines on the same axis
    Xs, and Ys, are lists of lists, corresponding to the same index in xlabels and ylabels
    data_labels is also a list corresponding to the Xs and Ys
    """
    # Define the font properties
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 12,
           }
    legend_font = {'family': 'serif', 'weight': 'normal', 'size': 8}
    plt.figure(figsize=figsize)  # Set the figure size
    if colors is None:
        for X, Y, label in zip(Xs, Ys, data_labels):
            plt.plot(X, Y, label=label, linewidth=linewidth) 
    else:
        for X, Y, label, color in zip(Xs, Ys, data_labels, colors):
            plt.plot(X, Y, label=label, linewidth=linewidth, color=color) 
    ax = plt.gca()
    # Set labels with the defined font properties
    plt.xlabel(xlabel, fontdict=font)
    plt.ylabel(ylabel, fontdict=font)
    # Set axis limits
    plt.xlim(xlim)
    plt.ylim(ylim)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Set tick label sizes
    plt.xticks(fontsize=font['size']-2)
    plt.yticks(fontsize=font['size']-2)
    ax.legend(prop=legend_font, frameon=False, bbox_to_anchor=legendloc, loc='upper left')
    ax.xaxis.set_major_locator(MultipleLocator(x_tick))
    # Save and close the plot
    plt.tight_layout()
    plt.savefig(f"{output_dir}{basename}_plt.png")
    plt.close()

def violin_plot(violin_df, output_file='/home/jacob.fenster/test.png', 
                x='experiment', y='score', hue='control', 
                marker_size=4):
    """ use sns here instead of matplolib https://stackoverflow.com/questions/62278350/matplotlib-seaborn-violin-plots-for-different-data-sizes
    Creates a violin plot and stripplot overlay
    violin_df is a df formatted with violin_master_formatter()
    x is 
    """ 
    fontsize = 18
    sns.set(style="whitegrid")
    plt.figure(figsize=(10,6))
    if hue == None:
        ax = sns.violinplot(data=violin_df, x=x, y=y, cut=0)
        ax = sns.stripplot(data=violin_df, x=x, y=y, jitter=True, dodge=True, color='black', size=marker_size, edgecolor='white', linewidth=0.2)
    else:
        ax = sns.violinplot(data=violin_df, x=x, y=y, hue=hue, cut=0)
        ax = sns.stripplot(data=violin_df, x=x, y=y, hue=hue, jitter=True, dodge=True, color='black', size=marker_size, edgecolor='white', linewidth=0.2)
    ax.set_xlabel(x, fontsize=fontsize)
    ax.set_ylabel(y, fontsize=fontsize)
    ax.grid(False)
    sns.despine()
    for spine in ax.spines.values():
        spine.set_color('black')
        spine.set_linewidth(1.5)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.savefig(output_file)

def stripplot(df, x='x', y='y', hue=None, output_dir='/home/jacob.fenster/', tag='exp_tag'):
    breakpoint()
    plt.figure(figsize=(10,6))
    if hue is None:
        ax = sns.stripplot(data=df, x=x, y=y, jitter=True, dodge=True, color='black', size=2, edgecolor='white', linewidth=0.2)
    else:
        ax = sns.stripplot(data=df, x=x, y=y, hue=hue, jitter=True, dodge=True, color='black', size=2, edgecolor='white', linewidth=0.2)
    plt.savefig(f"{output_dir}{tag}_plt.png")


def four_plot_fig_pmsa_summary(num_alns, avg_IDs, avg_covs, N_80, Neff, 
                               out_dir, basename):
    """
    Plots four figures summarizing the pMSAs of a given organism.
    Top left figure is hist of # alignments(num_alngs), 
    Top right is hist of avg %id (avg_IDs) and avg %coverage (avg_covs)
    Bottom left is hist of number of sequences after clustering at 80% id (N_80)
    bottom right is hist of Neff_80, N_80 / aln_length^0.5 (Neff)
    num_alns, avg_IDs, avg_covs, N_80, and Neff are lists
    Saves to workding dir with {out_dir}{basename}_runsummaryplt.png

    num_alns = summary['#_alns'].tolist()
    avg_IDs = summary['mean_ID'].tolist()
    avg_covs = summary['mean_cov'].tolist()
    N_80 = summary['Nseq_80'].tolist()
    Neff = summary['Neff_80'].tolist()
    """
    tick_fontsize, label_fontsize = 22, 24
    num_values = len(avg_IDs)
    if num_values < 100:
        n_bins = 3 * int(math.sqrt(num_values))
    else: 
        n_bins = int(math.sqrt(num_values))
        
    # Create a figure and a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 12))
    
    # Add an overall figure title
    title_basename = os.path.basename(basename).split('-pMSA')[0]
    fig.suptitle(f'{title_basename}', fontsize=label_fontsize+2)
    
    # Top-left subplot: Histogram of #_alns
    axs[0, 0].hist(num_alns, bins=n_bins, color='g', density=True)
    axs[0, 0].set_xlim([0, None])
    axs[0, 0].set_xlabel('#_alns', fontsize=label_fontsize)
    axs[0, 0].set_ylabel('Frequency', fontsize=label_fontsize)
    current_ticks = axs[0, 0].get_xticks()
    new_ticks = current_ticks[::2]  # take every second tick
    axs[0, 0].set_xticks(new_ticks)

    # Top-right subplot: Histogram of %ID and %Coverage
    axs[0, 1].hist(avg_IDs, bins=n_bins, alpha=0.5, label='average %ID', color='b', density=True)
    axs[0, 1].hist(avg_covs, bins=n_bins, alpha=0.6, label='average %Coverage', color='#F28C28', density=True)
    axs[0, 1].set_xlabel('avg %ID or avg %cov', fontsize=label_fontsize)
    axs[0, 1].set_xlim([0, 100])
    axs[0, 1].legend(loc='upper center', fontsize=(tick_fontsize-4), frameon=False)
    
    # Bottom-left subplot: Histogram of N_80
    axs[1, 0].hist(N_80, bins=n_bins, color='r', density=True)
    axs[1, 0].set_xlim([0, None])
    axs[1, 0].set_xlabel('N_80', fontsize=label_fontsize)
    axs[1, 0].set_ylabel('Frequency', fontsize=label_fontsize)
    
    # Bottom-right subplot: Histogram of Neff
    axs[1, 1].hist(Neff, bins=n_bins, color='m', density=True)
    axs[1, 1].set_xlim([0, None])
    axs[1, 1].set_xlabel('Neff', fontsize=label_fontsize)
    axs[1, 1].set_ylabel('Frequency', fontsize=label_fontsize)
    
    # Update axis tick fonts and style for all subplots
    for i in range(2):
        for j in range(2):
            for tick in axs[i, j].xaxis.get_major_ticks():
                tick.label1.set_fontsize(tick_fontsize)
            for tick in axs[i, j].yaxis.get_major_ticks():
                tick.label1.set_fontsize(tick_fontsize)
            axs[i, j].spines['top'].set_visible(False)
            axs[i, j].spines['right'].set_visible(False)
            axs[i, j].spines['left'].set_linewidth(2)
            axs[i, j].spines['bottom'].set_linewidth(2)
            axs[i, j].tick_params(axis='both', width=2, size=6)

    #plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.savefig(f"{out_dir}{basename}_runsummaryplt.png")
    plt.close()


def plot_scatter(X, Y, out_file,
                 x_label='x', y_label='y',
                 xlim=(0,1), ylim=(0,1), tick_interval=0.2, plot_size=(5,5), 
                 plot_type='ppt'):
    # Creates a scatter plot with clickable data points that will auto annotate
    # X and Y are lists of data. x_ y_label are axis label strings
    # xlim and ylim are tuples of length 2 to set axis limit
    # custom_labels is a list of data point labels whose index corresponds to X and Y
    plot_types = ['ppt']
    if plot_type not in plot_types:
        print(f"ERROR: plot type {plot_type} not an option. Exiting...")
        sys.exit()
    if plot_type == 'ppt':
        font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 18}
        tick_width, tick_length = 1, 5
        point_size, point_color = 3, 'black' 
    print(f"Plotting {plot_type} scatter with:\nxlabel={x_label}, ylabel={y_label}, fig_size={plot_size}\n\
xlim={xlim}, ylim={ylim}, tick interval={tick_interval}")
    plt.figure(figsize=plot_size)
    plt.scatter(X, Y, s=point_size, c=point_color, marker='.')
    plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
    plt.xlabel(x_label, fontdict=font)
    plt.ylabel(y_label, fontdict=font)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.tick_params(axis='both', labelsize=font['size'], direction='out', width=tick_width, length=tick_length, labelfontfamily=font['family'])
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()

def plot_scatter_linear_reg(X, Y, out_file,
                 x_label='x', y_label='y',
                 xlim=(0,1.05), ylim=(0,1.05), tick_interval=0.2, plot_size=(5,5), 
                 xscale=None, yscale=None, lin_reg=True, plot_type='ppt_bigdata'):
    # Creates a scatter plot with clickable data points that will auto annotate
    # X and Y are lists of data. x_ y_label are axis label strings
    # xlim and ylim are tuples of length 2 to set axis limit
    # custom_labels is a list of data point labels whose index corresponds to X and Y
    plot_types = ['ppt_bigdata', 'ppt_smalldata']
    if plot_type not in plot_types:
        print(f"ERROR: plot type {plot_type} not an option. Exiting...")
        sys.exit()
    if plot_type == 'ppt_bigdata' or plot_type == 'ppt_smalldata':
        font = {'family': 'serif',
                'color':  'black',
                'weight': 'normal',
                'size': 18}
        tick_width, tick_length = 1, 5
        point_size, point_color = 3, 'black' 
        if plot_type == 'ppt_smalldata':
            point_size = 12
    print(f"Plotting {plot_type} scatter with:\nxlabel={x_label}, ylabel={y_label}, fig_size={plot_size}\n\
xlim={xlim}, ylim={ylim}, tick interval={tick_interval}")
# linear regression
    # Clean out nan values
    if type(X) != np.ndarray and type(Y) != np.ndarray:
        X, Y = np.array(X), np.array(Y)
    mask = ~np.isnan(X) & ~np.isnan(Y)
    X_clean, Y_clean = X[mask], Y[mask]
    if lin_reg:
        # calculate linear best fit
        coefficients = np.polyfit(np.array(X_clean), np.array(Y_clean), 1, full=False)
        # calculate residual sum of squares and total sum of square
        slope, intercept = coefficients
        y_pred = np.polyval(coefficients, X_clean)
        residuals = Y_clean - y_pred
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((Y_clean-np.mean(Y_clean))**2)
        r_squared = 1 - (ss_res/ss_tot)
        # plot the fit line on the same scatter 
        line = slope * X_clean + intercept
    plt.figure(figsize=plot_size)
    ax = plt.gca()
    plt.scatter(X_clean, Y_clean, s=point_size, c=point_color, marker='.', label='data')
    plt.xlim(xlim)
    plt.ylim(ylim)
    if lin_reg:
        plt.plot(X_clean, line, linestyle='-', marker=None, color='red')
        plt.text(.98, 1.07, r'$R^2$' + f':{r_squared:.2f}', fontsize=font['size']-2, color='black', ha='right', va='top', transform=ax.transAxes)
    if tick_interval is not None:
        plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
        plt.gca().yaxis.set_major_locator(ticker.MultipleLocator(tick_interval))
        plt.tick_params(axis='both', labelsize=font['size']-2, direction='out', width=tick_width, length=tick_length, labelfontfamily=font['family'])
    if xscale == 'log':
        plt.xscale('log')
        plt.gca().xaxis.set_major_locator(LogLocator(base=10, numticks=10))
        plt.tick_params(axis='x', labelsize=font['size']-2, direction='out', width=tick_width, length=tick_length, labelfontfamily=font['family'])
        plt.tick_params(axis='y', labelsize=font['size']-2, direction='out', width=tick_width, length=tick_length, labelfontfamily=font['family'])
    if yscale == 'log':
        plt.yscale('log')
        plt.gca().yaxis.set_major_locator(LogLocator(base=10, numticks=10))
        plt.tick_params(axis='y', labelsize=font['size']-2, direction='out', width=tick_width, length=tick_length, labelfontfamily=font['family'])
    plt.xlabel(x_label, fontdict=font)
    plt.ylabel(y_label, fontdict=font)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(out_file)
    plt.close()



# Below are three functions that worked at one point to generate scatter plots that allowed for interactive 
# clicking of data points to label them. Not working as of 3/12/24
def _create_annot(ax, sc, labels, ind):
    # creates a clickable annotation for the interactive scatter
    pos = sc.get_offsets()[ind["ind"][0]]
    annot = ax.annotate("", xy=pos, xytext=(10,10), textcoords="offset points",
                        arrowprops=dict(arrowstyle="->", lw=0.5), fontsize=8)
    text = ', '.join([labels[n] for n in ind['ind']])
    annot.set_text(text)
    return annot

def _onclick(event, ax, sc, annot_list, labels):
    # creates a clickable event of a data point in a scatter plot 
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
            annot = _create_annot(ax, sc, labels, ind)
            annot.set_visible(True)
            annot_list.append(annot)
            plt.draw()

def plot_clickable_scatter(X, Y, custom_labels,  x_label, y_label, xlim, ylim):
    # Creates a scatter plot with clickable data points that will auto annotate
    # X and Y are lists of data. x_ y_label are axis label strings
    # xlim and ylim are tuples of length 2 to set axis limit
    # custom_labels is a list of data point labels whose index corresponds to X and Y
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
    fig.canvas.mpl_connect("button_press_event", lambda event: _onclick(event, ax, sc, annot_list, custom_labels))

    plt.show()

def plot_single_hist_density(df, col, output_file, 
                             xlim=(0,100), xlabel='xlabel', n_bins=None,
                             title=None):
    """
    Plots histograms for specififed column
    - df (DataFrame): The input DataFrame containing the data
    - column_names (list): List of column name for which to plot histogram
    """
    # Calculate statistics
    mean_value, median_value, stdev, num_values = df[col].mean(), df[col].median(), df[col].std(), df[col].count()
    if n_bins is None:
        if num_values < 100:
            n_bins = 3 * int(np.sqrt(len(df[col])))
        else: 
            n_bins = int(np.sqrt(len(df[col])))
    plt.figure(figsize=(6, 6))
    plt.hist(df[col], bins=n_bins, edgecolor='black', alpha=0.7, density=True)
    plt.xlim(xlim)
    plt.rcParams.update({'font.size': 20})
    plt.annotate(f"Mean: {mean_value:.1f}\nMedian: {median_value:.1f}\nStd Dev: {stdev:.1f}\nn: \
            {num_values}\nn_bins: {n_bins}", xy=(0.6, 0.8), xycoords='axes fraction', fontsize=12)
    if title is not None:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Relative Frequency')
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return

def plot_single_hist(df, col, output_file, 
                     xlim=(0,100), xlabel='xlabel', n_bins=None,
                     title=None, yaxis='linear'):
    """
    Plots histograms for specififed column
    - df (DataFrame): The input DataFrame containing the data
    - column_names (list): List of column name for which to plot histogram
    """
    # Calculate statistics
    mean_value, max_value, stdev, num_values = df[col].mean(), df[col].max(), df[col].std(), df[col].count()
    if n_bins is None:
        if num_values < 100:
            n_bins = 3 * int(np.sqrt(len(df[col])))
        else: 
            n_bins = int(np.sqrt(len(df[col])))
    plt.figure(figsize=(5, 4))
    plt.hist(df[col], bins=n_bins, edgecolor='black', alpha=0.7)
    if yaxis == 'log':
        plt.yscale('log')
    ax = plt.gca()
    plt.xlim(xlim)
    plt.rcParams.update({'font.size': 12})
    plt.annotate(f"Mean: {mean_value:.1f}\nMax: {max_value:.1f}\nn: {num_values}", 
                 xy=(0.6, 0.8), xycoords='axes fraction', fontsize=10)
    if title is not None:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Save the plot
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    return

def plot_pie_chart(data, labels, column='',
                  summary_text=(.01,.01), figsize=(3,3), output_file=''):
    """
    data is list of floats/integers
    lables is list of strings associated with data list
    """
    plt.figure(figsize=figsize)
    plt.pie(data, labels=labels, startangle=90)
    plt.axis('equal')
    plt_text = f'n {column} = {len(labels)}'
    plt.text(summary_text[0], summary_text[1], plt_text, fontsize=10)
    plt.savefig(output_file)
