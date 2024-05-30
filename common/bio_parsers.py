# This contains functions to parse various file types containing biological data
# last updated 20240119
import os
import pandas as pd
from typing import Dict, List, Tuple, Set

def read_fasta(file_path: str):
    """
    Reads a FASTA file and returns a list of descs and a list of seqs.
    Arguments:
        file_path: path to fasta formatted file
    Returns:
        two lists, descs and seqs lists which contain the desription and sequence 
        of each fasta entry with same indicies
    """
    descs = []
    seqs = []
    sequence = ""
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()  # Remove leading and trailing whitespaces
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):  # Description line
                if sequence:
                    seqs.append(sequence)
                    sequence = ""
                descs.append(line[1:])
            else:
                sequence += line
        if sequence:
            seqs.append(sequence)
    return descs, seqs

def read_hmmer_sto(input_sto):
#this reads a Stockholm 1.0 format MSA and parses between blocks and allows for redundant hits
# Uses '#=GC RF' as end of block identifier. This is likely HMMER specific 
# jackhmmer outputs the query alignemnt as the first line in the alignments
# if using phmmer set to true to remove leading hit counter ####|desc
    query, query_id = True, '' #holder to find query for block parsing
    unique = False # this is to check if there is the warning for unique seq IDs
    seq_ids = [] 
    #first extract the query seqeunce and id list of all alignments
    with open(input_sto, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('#') or line.startswith('\n'):
            if '# WARNING: seq names have been made unique by adding a prefix of \"<seq#>|\"' in line:
                unique = True
            elif line.startswith('#=GC RF'):
                break #stop after first block, should work if only one block
            continue
        else:
            seq_id, seq = line.strip().split(None, 1)
            if unique:
                unique_id = ""
                for chunk in seq_id.split('/')[0].split('|')[1:]:
                   unique_id = unique_id + '|' +  chunk
                seq_ids.append(unique_id[1:])
            else:
                seq_ids.append(seq_id.split('/')[0])
            if query: #record first entry as query
                query_id, query = seq_id, False
    seqs = ['' for _ in range(len(seq_ids))]
    for i in range(0, len(lines)):
        if lines[i].startswith(query_id): #parse block
            l = i #pull line index for use in while loop
            j = 0 #initialize hit alignment counter
            block = True #loop intil end of block line
            while block:
                if lines[l].startswith('#=GC RF'):
                    block = False
                if lines[l].startswith('#'):
                    l += 1 #move to next line
                    continue
                else:
                    seq_id, seq = lines[l].strip().split(None, 1)
                    seqs[j] += seq
                    l += 1
                    j += 1
    return seq_ids, seqs

def fasta_out_single_line(descs, seqs, output_path):
    # this outputs a fasta/a3m file output where sequences are output onto a single line regardless of length
    # adds an extra line after the last sequence
    with open(output_path, 'w') as f:
        for desc, seq in zip(descs, seqs):
            f.write(f">{desc}\n{seq}\n")

def parse_uniprot_proteome(proteome_file):
    # this parses the accession number for each entry in a UniProt proteome file and returns a list
    proteome_df = pd.DataFrame(columns=['description', 'sequence'])
    proteome_df.index.name = 'uniprot_accession'
    descs, seqs = read_fasta(proteome_file)
    for desc, seq in zip(descs, seqs):
       accession = desc.split('|')[1]
       proteome_df.loc[accession] = {'description':desc, 'sequence':seq}
    return proteome_df 

