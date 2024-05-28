### pMSA analysis functions to calculate metrics of pMSAs
### This script is to calculate the Neff80 values of pMSA inputs and load to dataframes 
import os, sys, glob, json, pdb
import pickle
import numpy as np
import pandas as pd 
# new to the AFppi 
import subprocess, math, statistics, time
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar


# GLOBAL variables
hhfilter = '/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter'
reformat = '/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl'

def _output_a3m(descs, seqs, basename, tmp):
    """
    outputs a temp fasta file with the input descs seqs
    file path is f"{tmp}{basename}"
    """
    with open(f"{tmp}{basename}", 'w') as f:
        for desc, seq in zip(descs, seqs):
            f.write(f">{desc}\n{seq}\n") 
    return f"{tmp}{basename}"

def _calculate_neff(pmsa,
                   id_tr=80, 
                   tmp='/home/jacob.fenster/tmp/',
                   hhfilter=hhfilter):
    """
    This calculates a custom Neff sequences using an hhfilter workaround. 
    The algorithm splits up the pMSA into the paired and unpaired sections, filters by the id_tr (default is 80%)
    splits the paired alignments by having a '__' motif.
    determines if the remaining alignments are unpaired protA or protB 
        by splitting by length and seeing if there are only gaps/alpha chars
    splits the two proteins by the amino acid lengths in the basename 'gene1-AA**__gene2-AA**-*.a3m'
    It outputs three Neff values for the pMSA, protA, and protB 
    
    pmsa is .a3m format pmsa file path. id_tr is percent id to filter at. default is 80%, 80
    tmp is the temp folder where the hhfilter file gets placed. This delets the filtered file after calculating
    hhfilter is path to the executible hhfilter binary
    """
    basename = os.path.basename(pmsa).split('.')[0]
    gene1 = basename.split('-AA')[0] 
    gene2 = basename.split('__')[1].split('-AA')[0]
    chain1 = int(basename.split('__')[0].split('-AA')[1])
    chain2 = int(basename.split('__')[1].split('-AA')[1].split('-')[0])
    descs, seqs = bpar.read_fasta(pmsa)
    paired_descs, paired_seqs = [], []
    chain1_descs, chain1_seqs = [], []
    chain2_descs, chain2_seqs = [], []
    for desc, seq in zip(descs, seqs):
        if all(char == '-' for char in seq[chain1:]):
            chain1_descs.append(desc), chain1_seqs.append(seq[:chain1])
        elif all(char == '-' for char in seq[:chain1]):
            chain2_descs.append(desc), chain2_seqs.append(seq[chain1:])
        elif '__' in desc: 
            paired_descs.append(desc), paired_seqs.append(seq)
        else:
            print(f"the desc: {desc} is strange.\n seq: {seq}")
    paired_msa = _output_a3m(paired_descs, paired_seqs, f'{basename}_paired.a3m', tmp)
    chain1_msa = _output_a3m(chain1_descs, chain1_seqs, f'{basename}_chain1.a3m', tmp)
    chain2_msa = _output_a3m(chain2_descs, chain2_seqs, f'{basename}_chain2.a3m', tmp)
    pmsa_command = f"{hhfilter} -v 0 -id {id_tr} -i {paired_msa} -o {tmp}{basename}-pfil{id_tr}.a3m"
    chain1_command = f"{hhfilter} -v 0 -id {id_tr} -i {chain1_msa} -o {tmp}{basename}-ch1fil{id_tr}.a3m"
    chain2_command = f"{hhfilter} -v 0 -id {id_tr} -i {chain2_msa} -o {tmp}{basename}-ch2fil{id_tr}.a3m"
    result_pmsa = subprocess.run(pmsa_command, shell=True, text=True, capture_output=True)
    result_chain1 = subprocess.run(chain1_command, shell=True, text=True, capture_output=True)
    result_chain2 = subprocess.run(chain2_command, shell=True, text=True, capture_output=True)

    if result_pmsa.returncode == 0:    
        descs, seqs = bpar.read_fasta(f"{tmp}{basename}-pfil{id_tr}.a3m")
        nseq_paired = len(descs)
    else:
        nseq_paired = 0

    if result_chain1.returncode == 0:
        descs, seqs = bpar.read_fasta(f"{tmp}{basename}-ch1fil{id_tr}.a3m")
        nseq_chain1 = len(descs)
    else:
        nseq_chain1 = 0

    if result_chain2.returncode == 0:
        descs, seq = bpar.read_fasta(f"{tmp}{basename}-ch2fil{id_tr}.a3m")
        nseq_chain2 = len(descs)
    else:
        nseq_chain2 = 0

    neff_paired = nseq_paired/math.sqrt(chain1+chain2)
    neff_chain1 = nseq_chain1/math.sqrt(chain1)
    neff_chain2 = nseq_chain2/math.sqrt(chain2)
    if gene1 == gene2: #return np.nan if homodimer set
        return np.nan, neff_chain1, neff_chain2
    else: 
        return neff_paired, neff_chain1, neff_chain2 


def _load_row_onto_master(values, 
                     columns, 
                     ppi,  
                     master):
    """
    Loads the value list into column list on the master dataframe given a ppi name
    searches for ppi match with columns 'PPI_A' and 'PPI_B'
    """
    found = False
    for index in master.index:
        if ppi == master.loc[index, 'PPI_A'] or ppi == master.loc[index, 'PPI_B']:
            found = True
            for value, column in zip(values, columns):
                master.loc[index, column] = value 
    if not found:
        print(f"{ppi} not found in master. Returning master")
        return master 
    else:
        return master

def load_neff(pmsa_dir,
              master,
              name_switch=None,
              id_tr=80, 
              tmp='/home/jacob.fenster/tmp/',
              hhfilter=hhfilter):
    """
    This calculates the Neff of a given pMSA and loads it onto the master df
    pmsa_dir is the directory of *.a3m format pMSAs
    master is the master dataframe
    name_switch is the switch for *.a3m parsing to load into master. name_switch options: 
        'None' has format gene1-AA**__gene2-AA***.*
        'fmdv' has format *768_gene1-AA**__*768_gene2-AA***.
    id_tr is the percent ID to calculate Neff at using the hhfilter workaround
    tmp is where the workaround hhfilter files are stored
    hhfilter is the path to the hhfilter binary
    """
    os.makedirs(tmp, exist_ok=True)
    pmsas = glob.glob(f"{pmsa_dir}*.a3m")
    for pmsa in pmsas:
        # parse PPI_A/PPI_B name
        basename = os.path.basename(pmsa).split('.')[0]
        if name_switch == None:
            gene1 = basename.split('-AA')[0] 
            gene2 = basename.split('__')[1].split('-AA')[0]
            ppi = f"{gene1}__{gene2}"
        elif name_switch == 'fmdv':
            gene1 = basename.split('__')[0].split('768_')[1].split('-AA')[0]
            gene2 = basename.split('__')[1].split('768_')[1].split('-AA')[0]
            ppi = f"{gene1}__{gene2}"
        neff_paired, neff_chain1, neff_chain2 = _calculate_neff(pmsa, id_tr, tmp, hhfilter)
        # assure that the genes are linked to the correct proteinA and proteinB
        neff_dict = {gene1: neff_chain1, gene2:neff_chain2}
        protA, protB = min([gene1, gene2]), max([gene1, gene2])
        master = _load_row_onto_master([neff_paired, neff_dict[protA], neff_dict[protB]], 
                                       [f'Neff{id_tr}-paired', f'Neff{id_tr}-protA', f'Neff{id_tr}-protB'], 
                                        ppi, master) 
    # cleanup temp!
    os.system(f'rm -r {tmp}')
    return master


######################## new
# optimizing functions and adding additional stats to add to master table

# adding ID_cov_gap_score calculations 
def _query_id_cov_gap_score(aln, query_aln):
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
    pid = 100*(identical/aligned)
    pcov = 100*(aligned/query_len)
    pgap = 100*(gaps/query_len)
    return pid, pcov, pgap 

def _msa_summary_stats(fasta_msa, query_aln):
    """
    Input a path to a fasta formatted MSA and the query_aln string
    outputs num_aln, mean_id, and mean_cov
    """
    descs, alns = bpar.read_fasta(fasta_msa)
    pids, pcovs = [], []
    if len(descs) == 0: # no aligmnets in msa condition 
        num_aln, mean_pid, mean_pcov = 0, np.nan, np.nan
        return num_aln, mean_pid, mean_pcov
    else:
        for aln in alns:
            pid, pcov, pgap = _query_id_cov_gap_score(aln, query_aln)
            pids.append(pid), pcovs.append(pcov)
        num_aln, mean_pid, mean_pcov = len(descs), statistics.mean(pids), statistics.mean(pcovs)
        return num_aln, mean_pid, mean_pcov
        
def _pmsa_stats(pmsa,
                name_switch=None, 
                id_tr=80, 
                tmp='/home/jacob.fenster/tmp/',
                hhfilter=hhfilter,
                reformat=reformat):
    """
    This calculates summary stats of a given pmsa. 
    outputs a row dictionary to be added onto a master df with PPI_A and PPI_B formatting

    This calculates a custom Neff sequences using an hhfilter workaround. 
    The algorithm splits up the pMSA into the paired and unpaired sections, filters by the id_tr (default is 80%)
    splits the paired alignments by having a '__' motif.
    determines if the remaining alignments are unpaired protA or protB 
        by splitting by length and seeing if there are only gaps/alpha chars
    splits the two proteins by the amino acid lengths in the basename 'gene1-AA**__gene2-AA**-*.a3m'
    It outputs three Neff values for the pMSA, protA, and protB 
    
    pmsa is .a3m format pmsa file path. id_tr is percent id to filter at. default is 80%, 80
    tmp is the temp folder where the hhfilter file gets placed. This delets the filtered file after calculating
    hhfilter is path to the executible hhfilter binary
    """
    basename = os.path.basename(pmsa).split('.')[0]
    chain1 = int(basename.split('__')[0].split('-AA')[1])
    chain2 = int(basename.split('__')[1].split('-AA')[1].split('-')[0])
    if name_switch == None:
        gene1 = basename.split('-AA')[0] 
        gene2 = basename.split('__')[1].split('-AA')[0]
    elif name_switch == 'fmdv':
        gene1 = basename.split('__')[0].split('768_')[1].split('-AA')[0]
        gene2 = basename.split('__')[1].split('768_')[1].split('-AA')[0]
    descs, seqs = bpar.read_fasta(pmsa)
    # Split into the paired and two unpaired MSAs
    paired_descs, paired_seqs = [], []
    chain1_descs, chain1_seqs = [], []
    chain2_descs, chain2_seqs = [], []
    for desc, seq in zip(descs, seqs):
        if all(char == '-' for char in seq[chain1:]):
            chain1_descs.append(desc), chain1_seqs.append(seq[:chain1])
        elif all(char == '-' for char in seq[:chain1]):
            chain2_descs.append(desc), chain2_seqs.append(seq[chain1:])
        elif '__' in desc: 
            paired_descs.append(desc), paired_seqs.append(seq)
        else:
            print(f"the desc: {desc} is strange.\n seq: {seq}")
    # output to a3m files for hhfiltering
    paired_msa = _output_a3m(paired_descs, paired_seqs, f'{basename}_paired.a3m', tmp)
    chain1_msa = _output_a3m(chain1_descs, chain1_seqs, f'{basename}_chain1.a3m', tmp)
    chain2_msa = _output_a3m(chain2_descs, chain2_seqs, f'{basename}_chain2.a3m', tmp)
    # hhfilter at id_tr and deal with no sequence error
    pmsa_command = f"{hhfilter} -v 0 -id {id_tr} -i {paired_msa} -o {tmp}{basename}-pfil{id_tr}.a3m"
    chain1_command = f"{hhfilter} -v 0 -id {id_tr} -i {chain1_msa} -o {tmp}{basename}-ch1fil{id_tr}.a3m"
    chain2_command = f"{hhfilter} -v 0 -id {id_tr} -i {chain2_msa} -o {tmp}{basename}-ch2fil{id_tr}.a3m"
    result_pmsa = subprocess.run(pmsa_command, shell=True, text=True, capture_output=True)
    result_chain1 = subprocess.run(chain1_command, shell=True, text=True, capture_output=True)
    result_chain2 = subprocess.run(chain2_command, shell=True, text=True, capture_output=True)
    if result_pmsa.returncode == 0:    
        descs, seqs = bpar.read_fasta(f"{tmp}{basename}-pfil{id_tr}.a3m")
        nseq_paired = len(descs)
    else:
        nseq_paired = 0
    if result_chain1.returncode == 0:
        descs, seqs = bpar.read_fasta(f"{tmp}{basename}-ch1fil{id_tr}.a3m")
        nseq_chain1 = len(descs)
    else:
        nseq_chain1 = 0
    if result_chain2.returncode == 0:
        descs, seq = bpar.read_fasta(f"{tmp}{basename}-ch2fil{id_tr}.a3m")
        nseq_chain2 = len(descs)
    else:
        nseq_chain2 = 0
    # calculate Neff at id_tr
    neff_paired = nseq_paired/math.sqrt(chain1+chain2)
    neff_chain1 = nseq_chain1/math.sqrt(chain1)
    neff_chain2 = nseq_chain2/math.sqrt(chain2)
    # deal with homodimers, return neff_paired as nan
    if gene1 == gene2: 
        neff_paired, nseq_paired = np.nan, np.nan
    # convert split alignments to fasta format for additional metric calculation
    reformat_pmsa_command = f"{reformat} a3m fas {tmp}{basename}_paired.a3m {tmp}{basename}_paired.fasta -v 0"
    reformat_chain1_command = f"{reformat} a3m fas {tmp}{basename}_chain1.a3m {tmp}{basename}_chain1.fasta -v 0"
    reformat_chain2_command = f"{reformat} a3m fas {tmp}{basename}_chain2.a3m {tmp}{basename}_chain2.fasta -v 0"
    result_pmsa_reformat = subprocess.run(reformat_pmsa_command, shell=True, text=True, capture_output=True)
    result_chain1_reformat = subprocess.run(reformat_chain1_command, shell=True, text=True, capture_output=True)
    result_chain2_reformat = subprocess.run(reformat_chain2_command, shell=True, text=True, capture_output=True)
    if os.path.exists(f"{tmp}{basename}_paired.fasta"):
        paired_query = paired_seqs[0]
        # handle homodimer case which has one paired query alignment
        if gene1 == gene2:
            paired_num_aln, paired_mean_pid, paired_mean_pcov = 1, np.nan, np.nan
        else:
            paired_num_aln, paired_mean_pid, paired_mean_pcov = _msa_summary_stats(f"{tmp}{basename}_paired.fasta", paired_query)
    else: 
         paired_num_aln, paired_mean_pid, paired_mean_pcov = 0, np.nan, np.nan
    if os.path.exists(f"{tmp}{basename}_chain1.fasta"):
        chain1_query = paired_seqs[0][:chain1]
        chain1_num_aln, chain1_mean_pid, chain1_mean_pcov = _msa_summary_stats(f"{tmp}{basename}_chain1.fasta", chain1_query)
    else: 
        chain1_num_aln, chain1_mean_pid, chain1_mean_pcov = 0, np.nan, np.nan
    if os.path.exists(f"{tmp}{basename}_chain2.fasta"):
        chain2_query = paired_seqs[0][chain1:]
        chain2_num_aln, chain2_mean_pid, chain2_mean_pcov = _msa_summary_stats(f"{tmp}{basename}_chain2.fasta", chain2_query)
    else: 
        chain2_num_aln, chain2_mean_pid, chain2_mean_pcov = 0, np.nan, np.nan
    protA, protB = min([gene1, gene2]), max([gene1, gene2])
    ppi = f"{protA}__{protB}"
    data = {
            gene1: {f'Nseq{id_tr}': nseq_chain1, f'Neff{id_tr}': neff_chain1, 'Nseq': chain1_num_aln,
                    'mean_pid': chain1_mean_pid, 'mean_pcov': chain1_mean_pcov}, 
            gene2: {f'Nseq{id_tr}': nseq_chain2, f'Neff{id_tr}': neff_chain2, 'Nseq': chain2_num_aln,
                    'mean_pid': chain2_mean_pid, 'mean_pcov': chain2_mean_pcov}
           }
    row = {
            f'Nseq{id_tr}-paired': nseq_paired, f'Neff{id_tr}-paired': neff_paired, 'Nseq-paired': paired_num_aln,
            'mean_pid-paired': paired_mean_pid, 'mean_pcov-paired': paired_mean_pcov,
            f'Nseq{id_tr}-protA': data[protA][f'Nseq{id_tr}'], f'Neff{id_tr}-protA': data[protA][f'Neff{id_tr}'], 'Nseq-protA': data[protA]['Nseq'],
            'mean_pid-protA': data[protA]['mean_pid'], 'mean_pcov-protA': data[protA]['mean_pcov'],
            f'Nseq{id_tr}-protB': data[protB][f'Nseq{id_tr}'], f'Neff{id_tr}-protB': data[protB][f'Neff{id_tr}'], 'Nseq-protB': data[protB]['Nseq'],
            'mean_pid-protB': data[protB]['mean_pid'], 'mean_pcov-protB': data[protB]['mean_pcov']
           }
    return row, ppi

def pmsa_dir_summary_stats(pmsa_dir,
                           master,
                           name_switch=None,
                           id_tr=80, 
                           tmp='/home/jacob.fenster/tmp/', 
                           hhfilter=hhfilter,
                           reformat=reformat):
    """
    This function takes a directory to a pMSA and a master dataframe, calculates stats on each pMSA
    master MUST HAVE INDICIES OF PPI_A
    name_switch is used to parse the protA and protB from the pMSA basenames. Can be 'None' or 'fmdv'. 
    """
    if type(master.index[0]) != str:
        print(f"ERROR: incorrect format for master index for pmsa_dir_summary_stats(). Convert index to PPA_A.\nexiting...")
        sys.exit()
    os.makedirs(tmp, exist_ok=True)
    pmsas = glob.glob(f"{pmsa_dir}*.a3m")
    for pmsa in pmsas:
        row, ppi =  _pmsa_stats(pmsa, name_switch=name_switch, id_tr=id_tr, tmp=tmp, hhfilter=hhfilter, reformat=reformat)
        for column in row.keys():
            master.loc[ppi, column] = row[column]
    # clean up tmp
    os.system(f"rm -r {tmp}")
    return master


