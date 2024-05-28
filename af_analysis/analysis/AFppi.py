import os, sys, glob, json, pdb
import pickle
import numpy as np
import pandas as pd 
# These definitions are to help with data processing of alphafold PPI models from the 
# Ian2021 algorithm and the AF-multimer algorithm. 
def max_contact_terimini_correction(arr):
# this definition inputs a RF2t numpy array of contact probabilities and excludes the last ten amino acids
# of the first protein and the first ten amino acids of the second protein when calculating the max contact PPI
    arr[-10:, :] = 0 # set the last 10 residues on C terminus of protein 1 to zero
    arr[:, :10] = 0 # set the first 10 residues on N terminus of protein 2 to zero
    max_contact = np.max(arr)
    return max_contact

def dir_max_contact(input_dir, protein_len_df):
    """
    This inputs a directory that contains the output .npz files from Ian2021 and 
    a dataframe containing the protein name as index with a length column as 'length'
    extracts the max contact probability from each .npz file using the max_contact_terimini_correction
    which igonores the last 10 C residues and first 10 N residues
    returns an output dataframe with index name 'PPI' and data labeled 'max_contact_prob'
    """
    npz_files = glob.glob(f"{input_dir}*.npz")
    output_df = pd.DataFrame(columns=['max_contact_prob'])
    output_df.index.name = 'PPI'
    incorrect_shape = 0
    for npz in npz_files:
        basename = os.path.basename(npz).split('.')[0]
        index = basename.split('_info')[0]
        protein1, protein2 = basename.split('__')[0], basename.split('__')[1].split('_info')[0]
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

def AFmult_data_extract(exp_dir):
    """This extracts AF multimer data from the exp_dir sub directories of the parent input_dir
    returns a df with the PPIs as indices and data columns 'iptm+ptm', 'ptm', and 'iptm' 
    """
    exp_dirs = glob.glob(f"{exp_dir}*/")
    df = pd.DataFrame(columns=['iptm+ptm', 'iptm', 'ptm'])
    df.index.name="PPI"
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
    """
    input the positive and negative control  dataframe with data in 'column', a threshold float and this 
    will calculate the precision and recall at that threshold.
    the ratio of positive to negative control rows/values  must be determined before this algorithm as this gives raw values
    """
    pos_above_tr = len(pos_ctrl_df[pos_ctrl_df[f"{column}"] > threshold])
    neg_above_tr = len(neg_ctrl_df[neg_ctrl_df[f"{column}"] > threshold]) 
    total_above_tr = pos_above_tr + neg_above_tr
    total_positive = len(pos_ctrl_df.index)
    if total_above_tr == 0:
        precision, recall = np.nan, np.nan
    else:
        precision = pos_above_tr / total_above_tr
        recall = pos_above_tr / total_positive 
    return precision, recall

def precision_recall_curve(pos_df, neg_df, threshold_arr, column):
    """
    input the positive control and negative control df, the threshold array, and the column name of the df to analyze
    the ratio of positive to negative control rows/values  must be determined before this algorithm as this does not normalize the ratio of
    signal to noise
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

def custom_merge_any(df1, df2, index_col, data_col1, data_col2):
    # this def does a custom merge allowing for a union of two dataframes given the index_col column
    # the same column name for index_col (which is not the index but a column used as unique index)
    # data_col1 and data_col2 are lists of columns to include from each respective dataframe in the merged df. 
    if type(df1.index.tolist()[0]) != int:
        df1 = df1.reset_index()
    if type(df2.index.tolist()[0]) != int:
        df2 = df2.reset_index()
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

def parse_pdb_to_df(pdb_file_path):
    # Define the columns
    columns = ['atom_number', 'atom_name', 'residue_name', 'chain_id',
               'residue_number', 'x', 'y', 'z', 'occupancy', 'plddt', 'element']

    # Initialize a list to hold all atom data
    atom_data = []

    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                # Split the line by whitespace
                parts = line.split()
                if len(parts) < 12 and len(parts[4]) > 2:
                    if  parts[4][0].isalpha() and parts[4][1].isdigit(): # OF error with residues >1000                    
                        atom_data.append((
                            int(parts[1]), parts[2], parts[3], parts[4][0], int(parts[4][1:]),
                            float(parts[5]), float(parts[6]), float(parts[7]),
                            float(parts[8]), float(parts[9]), parts[-1]
                        ))
                    else:
                        breakpoint()
                        continue
                else:
                # Extract the required elements from the parts list
                # Adjust the indices if necessary based on the structure of your file
                    try:
                        atom_data.append((
                            int(parts[1]), parts[2], parts[3], parts[4], int(parts[5]),
                            float(parts[6]), float(parts[7]), float(parts[8]),
                            float(parts[9]), float(parts[10]), parts[-1]
                        ))
                    except ValueError:
                        print(f"incorrect format for {pdb_file_path}. skipping...")
                        return None
    # Create a DataFrame from the atom data
    df = pd.DataFrame(atom_data, columns=columns)
    return df

def avg_plddt_from_pdb(pdb_file_path):
    # input pdb file path that has the pDLLT confidence value in place of the b_factors and this parses with parse_pdb_to_df
    # and calculates and returns the average pDLLT value
    df = parse_pdb_to_df(pdb_file_path)
    if df is not None:
        avg_plddt = df['plddt'].mean()
    else:
        avg_plddt = np.nan
    return avg_plddt

def _return_AF_uniref90(protein, af_output_dir):
    # this inputs the parent directory of alphafold multimer output and finds the path of the uniref90 file
    # Returns None if didn't find the protein folder or the uniref90_hits.sto file
    MSA_dirs = glob.glob(f"{af_output_dir}*/")
    for dir in MSA_dirs:
        with open(f"{dir}msas/chain_id_map.json", 'r') as f:
            data = json.load(f)
        protA, protB = data['A']['description'].split('__')[0], data['B']['description'].split('__')[0]
        print(f"{protA}, {protB}")
        if protA == protein:
            if os.path.exists(f"{dir}msas/A/uniref90_hits.sto"):
                return f"{dir}msas/A/uniref90_hits.sto" 
            else:
                return None
        elif protB == protein:
            if os.path.exists(f"{dir}msas/B/uniref90_hits.sto"):
                return f"{dir}msas/B/uniref90_hits.sto" 
    return None

def afmult_template_search(protein_name, af_data_parentdir, temp='home/jacob.fenster/temp', reformat='/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl', hhsearch='hhsearch', pdb70='/lustrefs/public/databases/alphafold/pdb70/'):
    # This script takes the AF multimer parent dir and searches subdirectories for the protein_name. it then 
    # 
    os.makedirs(temp, exist_ok=True)
    uniref90_sto = _return_AF_uniref90(protein_name, af_data_parentdir)
    if uniref90_sto is None:
        print(f"No AFmult uniref90 sto file for {protein_name}. skipping...")
        return None
    else:
        # convert to a3m and conduct an hhsearch 
        os.system(f"{reformat} sto a3m {uniref90_sto} {temp}uniref90_hits-{protein_name}.a3m")
        os.system(f"{hhsearch} -i {temp}uniref90_hits-{protein_name}.a3m -d {pdb70}pdb70 -o {temp}{protein_name}-pdb70hits.hhr -maxres 1000000")
        # read hhsearch .hhr file to string
        with open(f"{temp}{protein_name}-pdb70hits.hhr", 'r') as f:
            hrr_string = f.read()
        return hrr_string

def _load_master_index(data_df, master, data_cols, master_cols):
    # Assumes matching indicies (PPI_A)
    # loads data_cols data from data_df onto master_cols master df
    # data_cols and master_cols are lists of column names
    for data_index in data_df.index:
        for data_column, master_column in zip(data_cols, master_cols):
            master.loc[data_index, master_column] = data_df.loc[data_index, data_column]
    return master
