# This package as functions that clean up and generate pMSAs from MSAs
# includes jackhmmer local alignment cleanup functions which handle duplicate hits per genome and protein
# and functions which build seed hmms for aligning sequences 
# filters sequences based on low coverage to query for output and \
# pairs sequences based on genome accession number for the NCBI all virus database.

import sys, os, glob, json, subprocess, pdb
import numpy as np
import pandas as pd
from collections import defaultdict
import common.bio_parsers as bpar 

def find_alifrom_alito(seq):
    alnlen = len(seq)
    alifrom, alito = None, None #index of start of alignment and end of alignment. First index of alignment is zero
    for c in range(alnlen): #iterate through alignment columns
        if seq[c].isalpha():
            if alifrom is None:
                alifrom = c
            alito = c #pulls the last alpha character in the list
    return alifrom, alito

def gap_ratio_score(seq, query_len):
#scores the gap ratio of seq versus query seq. 
# input is sequence and length of the query. So for gap removed sequences is the gap removed seq and len of query no gaps
    gap = 0
    for c in range(len(seq)):
        if seq[c] == '-':
            gap += 1
    gap_ratio = gap/query_len
    return gap_ratio

def ID_score(seq, query_seq):
    aligned, identical = 0, 0
    for c in range(len(query_seq)):
        if seq[c].isupper():
            aligned += 1
            if seq[c] == query_seq[c]:
                identical += 1
    ID = identical/aligned
    return ID

def overlapping_ID_score(seqin, overlap, query_seqin):
    #scores the percent ID of the given overlapping region. Input is full alignment of seq and query. this slices overlap
    seq, query_seq = seqin[overlap[0]:overlap[1]+1], query_seqin[overlap[0]:overlap[1]+1]
    identical, aligned = 0, 0
    for c in range(len(query_seq)):
        if seq[c].isupper():
            aligned +=1
            if seq[c] == query_seq[c]:
                identical += 1
    try:
        ID = identical/aligned
    except ZeroDivisionError:
        return 0
    return ID

def coverage_score(seq, query_len):
    aligned = 0
    for c in range(len(seq)):
        if seq[c].isupper():
            aligned += 1
    cov = aligned/query_len
    return cov

def find_duplicates(seq_ids, seqs):
    #record indicies in original list
    entry_to_indices = defaultdict(list)
    for i, entry in enumerate(seq_ids):
        entry_to_indices[entry].append(i)
    # Find duplicates
    duplicates = {entry: indices for entry, indices in entry_to_indices.items() if len(indices) > 1}
    #Find alignment to and from for each duplicated alignment and add to dictionary of duplicates (alito, alifrom, seq)
    for key in duplicates:
        for i in range(len(duplicates[key])):
            alifrom, alito = find_alifrom_alito(seqs[duplicates[key][i]])
            duplicates[key][i] = (alifrom, alito, seqs[duplicates[key][i]]) #duplicates[key][i] = (alifrom, alito, sequence)
    return duplicates

def remove_duplicate_descs(descs, seqs):
    # this inputs a list of descs and sequences and removes any duplicated names, taking the first occurance of each name
    # also removes any unknown amino acids
    # intended to be used for concatinating af data 
    descs_out, seqs_out = [], []
    unique = {}
    for desc, seq in zip(descs, seqs):
        if desc in unique:
            continue
        elif any(char in seq for char in "xXbBjJzZ"):
            continue
        else:
            unique[desc] = 1
            descs_out.append(desc), seqs_out.append(seq)
    return descs_out, seqs_out

def remove_duplicates_prioitize_NCBI(seq_ids, seqs):
    # this one inputs a list of descs and sequences and finds idential sequences and 
    # creates a new list with only one exmaple of each unique sequence
    # prioritizes sequences from the NCBI all virus database ie. 'lcl|' header 
    unique_seq_ids, unique_seqs = [], []
    query_id = seq_ids[0]
    duplicates = {} 
    # initialize the dictionary 
    # dictionary has strucutre {unique alignment:[list of descriptions], resolved? T/F]}
    for seq in seqs:
        duplicates[seq] = [[],False]
    # load dictionary
    for seq_id, seq in zip(seq_ids, seqs):
        duplicates[seq][0] = duplicates[seq][0] + [seq_id]
    # generate output keeping only one of each duplicate prioritizing NCBI sequences
    seq_ids_out, seqs_out = [], [] 
    for seq_id, seq in zip(seq_ids, seqs):
        if len(duplicates[seq][0]) > 1:
            for duplicate_id in duplicates[seq][0]:
                if (duplicate_id.split('|')[0] == 'lcl' or duplicate_id == query_id) and not duplicates[seq][1]: #this will include multiples from 'lcl' if they exist.
                    seq_ids_out.append(duplicate_id), seqs_out.append(seq)
            for duplicate_id in duplicates[seq][0]:
                if (duplicate_id.split('|')[0] == 'lcl' or duplicate_id == query_id) and not duplicates[seq][1]: #this will include multiples from 'lcl' if they exist.
                    duplicates[seq][1] = True
            if not duplicates[seq][1]:
                # pick first description if no 'lcl'
                seq_ids_out.append(duplicates[seq][0][0]), seqs_out.append(seq)
                duplicates[seq][1] = True
        else:
            seq_ids_out.append(seq_id), seqs_out.append(seq)
    return seq_ids_out, seqs_out
                    
def resolve_pairwise_overlaps(duplicates, query_seq):
# this inputs duplicate local MSAs, finds the type of overlap between all pairwise combinations
# and concatinates sequences, keeping the highest %ID if overlapping. Returns same dictionary with 
# one round of consolidation.  
    for key in duplicates:
        pairs = []
        #pair adjacent and leave lone last entry
        for i in range(0,len(duplicates[key]), 2):
            try:
                pairs.append((duplicates[key][i], duplicates[key][i+1]))
            except IndexError:
                pairs.append((duplicates[key][i], duplicates[key][i])) #filler for odd out
        duplicates[key] = []
        for pair in pairs:
            concat_seq = ''
            #test if overlapping and make overlap tuple (start, finish, index1, index2) or (None, None, i1, i2) if no overlap
            if (pair[0][0]>pair[1][0]) and (pair[0][1]<pair[1][1]): #Case1, 1 contained in 2 
                overlap = (pair[0][0], pair[0][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[1][2][0:pair[0][0]] + pair[0][2][pair[0][0]:pair[0][1]+1] + pair[1][2][pair[0][1]+1:]
                else:
                    concat_seq = pair[1][2]
            elif (pair[0][0] < pair[1][0]) and (pair[0][1] > pair[1][1]): #Case2, 2 contained in 1 
                overlap = (pair[1][0], pair[1][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[0][2]
                else:
                    concat_seq = pair[0][2][0:pair[1][0]] + pair[1][2][pair[1][0]:pair[1][1]+1] + pair[0][2][pair[1][1]+1:]
            elif (pair[0][0]<pair[1][0]) and (pair[0][1]>pair[1][0]) and (pair[0][1]<pair[1][1]): #Case3, end of 1 overlaps with start of 2 *!
                overlap = (pair[1][0], pair[0][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[0][2][0:pair[0][1]+1] + pair[1][2][pair[0][1]+1:]
                else:
                    concat_seq = pair[0][2][0:pair[1][0]] + pair[1][2][pair[1][0]:]
            elif (pair[0][0]>pair[1][0]) and (pair[0][0]<pair[1][1]) and (pair[0][1]>pair[1][1]): #Case4, end of 2 overlaps with start of 1
                overlap = (pair[0][0], pair[1][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[1][2][0:pair[0][0]] + pair[0][2][pair[0][0]:] 
                else:
                    concat_seq = pair[1][2][0:pair[1][1]+1] + pair[0][2][pair[1][1]+1:]
            elif (pair[0][0]==pair[1][0]) and (pair[0][1]<pair[1][1]): #Case 5, if same start but 2 is longer
                overlap = (pair[0][0], pair[0][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[0][2][0:pair[0][1]+1] + pair[1][2][pair[0][1]+1:]
                else:
                    concat_seq = pair[1][2]
            elif (pair[0][0]==pair[1][0]) and (pair[0][1]>pair[1][1]): #Case 6, same start but 1 is longer
                overlap = (pair[1][0], pair[1][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[0][2]
                else:
                    concat_seq = pair[1][2][0:pair[1][1]+1] + pair[0][2][pair[1][1]+1:]
            elif (pair[0][0] > pair[1][0]) and (pair[0][1]==pair[1][1]): #Case7, if same finish but 2 is longer
                overlap = (pair[0][0], pair[0][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[1][2][0:pair[0][0]] + pair[0][2][pair[0][0]:]
                else:
                    concat_seq = pair[1][2]
            elif (pair[0][0] < pair[1][0]) and (pair[0][1]==pair[1][1]): #Case8, if same finish but 1 is longer
                overlap = (pair[1][0], pair[1][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1 > score2:
                    concat_seq = pair[0][2]
                else:
                    concat_seq = pair[0][2][0:pair[1][0]] + pair[1][2][pair[1][0]:]
            elif (pair[0][0]==pair[1][0]) and (pair[0][1]==pair[1][1]): #Case9, same start and finish
                overlap = (pair[0][0], pair[0][1])
                score1, score2 = overlapping_ID_score(pair[0][2], overlap, query_seq), overlapping_ID_score(pair[1][2], overlap, query_seq)
                if score1>score2:
                    concat_seq = pair[0][2]
                else:
                    concat_seq = pair[1][2]
            else: #if not overlapping
                if pair[0][0] < pair[1][0]:
                    #these are inclusive of the to, or finish character 
                    concat_seq = pair[0][2][0:pair[0][1]+1] + pair[1][2][pair[0][1]+1:]
                else: 
                    concat_seq = pair[1][2][0:pair[1][1]+1] + pair[0][2][pair[1][1]+1:]
            alito, alifrom = find_alifrom_alito(concat_seq)
            duplicates[key].append((alito, alifrom, concat_seq))
    return duplicates

def merge_all_overlaps(duplicates, query_seq):
    resolved = {key: None for key in duplicates.keys()}
    while len(duplicates) > 0:
        for key in resolved:
            try:
                if len(duplicates[key]) == 1:
                    resolved[key] = duplicates[key][0][2]
                    resolved[key] = [resolved[key], False] #add tag for not loaded back into MSA
                    del duplicates[key]
            except KeyError:
                continue
        duplicates = resolve_pairwise_overlaps(duplicates, query_seq)
    return resolved

def count_fasta(fasta_file):
#intputs a path to a fasta file and outputs the number of entries in the file 
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count

def phmmer_stats_resolved_duplicates(seq_ids, seqs, query_seq, resolved):
#This inputs the list of seq_IDs, the sequences, and the resolved list of duplicate hits from merge_all_overlaps
# It builds a dataframe only adding the resolved sequence for duplicates and calculates the desired stats on the sequences in the df
    query_len = float(len(query_seq.replace('.','').replace('-','')))
    df = pd.DataFrame(columns=['Hits', 'alignments', '%ID', 'coverage', 'gap_ratio'])
    for i in range(0, len(seq_ids)): #iterate through names
        if (seq_ids[i] in resolved) and (resolved[seq_ids[i]][1] == True):
            continue # skip this entry if duplicate has been added to dataframe already
        if (seq_ids[i] in resolved) and (resolved[seq_ids[i]][1] == False):
            seqs[i] = resolved[seq_ids[i]][0]
            resolved[seq_ids[i]][1] = True
        df.loc[i] = [seq_ids[i], seqs[i], ID_score(seqs[i], query_seq), coverage_score(seqs[i], query_len), gap_ratio_score(seqs[i], query_len)]
    return df

def phmmer_sim_filter(df, ID_tr, cov_tr, gap_ratio_tr):
# pulls rows that have values equal to or gretor than %ID and coverage thresholds, and less than gap ratio threshold
    filtered_df = df[(df['%ID']>ID_tr) & (df['coverage']>cov_tr) & (df['gap_ratio']<gap_ratio_tr)]
    rows = len(filtered_df.axes[0])
    return filtered_df, rows

def seed_sto_out(df, file_output):
#this inputs a pd DataFrame with 'Hits' column, 'alisngments' column and saves a STOCKHOLM 1.0 format file
    with open(file_output, 'w') as handle:
        handle.write("# STOCKHOLM 1.0\n")
        max_hit_length = df['Hits'].apply(len).max()
        padding = 4
        for index, row in df.iterrows():
            hit = row['Hits']
            alignment = row['alignments']
            # Calculate the number of spaces needed to align this entry
            num_spaces = max_hit_length - len(hit) + padding
            handle.write(f"{hit}{' ' * num_spaces}{alignment}\n")
        # Write the Stockholm footer
        handle.write("//\n")

def merge_overlaps_and_update(duplicates, query_seq, seq_ids, seqs):
    resolved = {key: None for key in duplicates.keys()}
    while len(duplicates) > 0:
        for key in resolved:
            try:
                if len(duplicates[key]) == 1:
                    resolved[key] = duplicates[key][0][2]
                    resolved[key] = [resolved[key], False] #add tag for not loaded back into MSA
                    del duplicates[key]
            except KeyError:
                continue
        duplicates = resolve_pairwise_overlaps(duplicates, query_seq)
    #update original lists of seqs and seq_ids.
    seq_ids_merge, seqs_merge = [], []
    for i in range(0, len(seq_ids)):
        if seq_ids[i] in resolved:
            if not resolved[seq_ids[i]][1]:
                resolved[seq_ids[i]][1] = True
                seq_ids_merge.append(seq_ids[i])
                seqs_merge.append(resolved[seq_ids[i]][0])
        else:
            seq_ids_merge.append(seq_ids[i])
            seqs_merge.append(seqs[i])
    return seq_ids_merge, seqs_merge

def highest_score_per_genome_update(seq_ids, seqs, query_seq):
    #this finds hits that are of different proteins within the same genome and picks the hit with the highest %ID score
    #ignores the query, first line, and keeps it
    genomes = []
    query = True
    for seq_id in seq_ids:
        if query:
            genomes.append(seq_id)
            query = False
            continue
        genome = seq_id.split('|')[1].split('.')[0]
        genomes.append(genome)
    #create dictionary of genomes with indicies in the original seq_ids and seqs
    entry_to_indices = defaultdict(list)
    for i, genome in enumerate(genomes):
        entry_to_indices[genome].append(i)
    #pull out dictionary of duplicated genome entries
    duplicates = {genome: indices for genome, indices in entry_to_indices.items() if len(indices) > 1}
    for genome in duplicates:
        ID_scores = []
        for i in duplicates[genome]: #iterate through indexes of duplicate genome hits
            score = ID_score(seqs[i], query_seq)
            ID_scores.append(score)
        index_of_max = ID_scores.index(max(ID_scores)) #this will return the first instance of the max score if there are duplicates
        winner_index = duplicates[genome][index_of_max]
        duplicates[genome] = [winner_index] #replace list of indicies with just the winnner
    seq_ids_merge, seqs_merge = [], []
    query = True
    for i in range(0, len(seq_ids)):
        if query:
            seq_ids_merge.append(seq_ids[i])
            seqs_merge.append(seqs[i])
            query = False
            continue
        genome = seq_ids[i].split('|')[1].split('.')[0]
        if genome in duplicates:
            if duplicates[genome][0] == i: #check if current index is equal to the max index of that genome
                seq_ids_merge.append(seq_ids[i])
                seqs_merge.append(seqs[i])
            else:
                continue
        else:
            seq_ids_merge.append(seq_ids[i])
            seqs_merge.append(seqs[i])
    return seq_ids_merge, seqs_merge


def gapfilter(seq_ids, seqs, query_seq,  gap_ratio_tr):
    # this filters out columns that are gaps in the query and removed alignments that are above the gap ratio tr after 
    # removing gapped columns 
    gap_indicies = [i for i, char in enumerate(query_seq) if char in ['-', '.']] #remove all - columns 
    gap_rem_query = ''.join([char for i, char in enumerate(query_seq) if i not in gap_indicies])
    gap_rem_descs, gap_rem_seqs = [], []
    for seq_id, seq in zip(seq_ids, seqs):
        gap_rem_seq = ''.join([char for i, char in enumerate(seq) if i not in gap_indicies])
        gap_rem_seq = gap_rem_seq.replace('.', '-') #replace insertion columns '.' by gap columns
        if gap_ratio_score(gap_rem_seq, len(gap_rem_query)) < gap_ratio_tr:
            gap_rem_descs.append(seq_id), gap_rem_seqs.append(gap_rem_seq)
    return gap_rem_descs, gap_rem_seqs

def gap_remove_and_filter(seq_ids, seqs, query_seq, resolved, gap_ratio_tr=0.5):
    # this takes in aligned seq_ids and seqs, the resolved dictionary of overlapped alignments
    # removes gap columns that are present in the query_seq, and filters those above the gap_ratio_tr
    seq_ids_out, seqs_out = [], []
    gap_indicies = [i for i, char in enumerate(query_seq) if char in ['-', '.']] #remove all - columns 
    gap_rem_query = ''.join([char for i, char in enumerate(query_seq) if i not in gap_indicies])
    for seq_id, seq in zip(seq_ids, seqs):
        if seq_id in resolved:
            if not resolved[seq_id][1]:
                gap_rem_seq = ''.join([char for i, char in enumerate(resolved[seq_id][0]) if i not in gap_indicies])
                gap_rem_seq = gap_rem_seq.replace('.', '-') #replace insertion columns '.' by gap columns
                if gap_ratio_score(gap_rem_seq, len(gap_rem_query)) < gap_ratio_tr:
                    seq_ids_out.append(seq_id), seqs_out.append(gap_rem_seq)
                resolved[seq_id][1] = True 
            continue
        gap_rem_seq = ''.join([char for i, char in enumerate(seq) if i not in gap_indicies])
        if gap_ratio_score(gap_rem_seq, len(gap_rem_query)) < gap_ratio_tr:
            seq_ids_out.append(seq_id), seqs_out.append(gap_rem_seq)
    return seq_ids_out, seqs_out

def sys_jackhmmer_search(fasta_query, seq_db, j_out, temp='/home/jacob.fenster/temp/', threads=1, iterations=3, jhmmer="hhsuite/jackhmmer"):
    # this def conducts a jackhmmer search using os.system() command 
    # returns the file path to the search result 
    # fasta_query - single protein query in fasta format. seq_db - fasta format database for search. j_out - where _jhmmer.sto files are saved
    os.makedirs(temp, exist_ok=True)
    basename = os.path.basename(fasta_query).split('.')[0]
    jhmmer_sto = f"{j_out}{basename}_jhmmer.sto"
    os.system(f'{jhmmer} -N {iterations} --cpu {threads} --notextw -o {temp}{basename}-jhmmer.log -A {jhmmer_sto} {fasta_query} {seq_db}')
    return jhmmer_sto

def sub_jackhmmer_search(fasta_query, seq_db, j_out, temp='/home/jacob.fenster/temp/', threads=1, iterations=3, jhmmer="hhsuite/jackhmmer"):
    # this def conducts a jackhmmer search using os.system() command 
    # returns the file path to the search result 
    # fasta_query - single protein query in fasta format. seq_db - fasta format database for search. j_out - where _jhmmer.sto files are saved
    os.makedirs(temp, exist_ok=True), os.makedirs(j_out, exist_ok=True)
    basename = os.path.basename(fasta_query).split('.')[0]
    jhmmer_sto = f"{j_out}{basename}_jhmmer.sto"
    command = f'{jhmmer} -N {iterations} --cpu {threads} --notextw -o {temp}{basename}-jhmmer.log -A {jhmmer_sto} {fasta_query} {seq_db}'
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode != 0:
        print(f"ERROR: jackhmmer search with {fasta_query} returncode {result.returncode}. stderr: {result.stderr}")
    else:
        print(f"jackhmmer search with {fasta_query} complete.")
        print(f"{result.stdout}")
        print("")
    return jhmmer_sto

def cleanup_sto_or_fasta(input_file, gap_ratio_tr=0.5, genome=False):
    # this inputs the path to a .sto alignment and cleans up for multiple hits per protein and genome
    # set genome=True to have this look for multiple hits per genome (all viral NCBI database jhmmer inputs)
    if os.path.basename(input_file).split('.')[1] == 'sto':
        descs, seqs = bpar.read_hmmer_sto(input_file)
    elif os.path.basename(input_file).split('.')[1] == 'fasta':
        descs, seqs = bpar.read_fasta(input_file)
    else:
        print(f"Incorrect input file")
        breakpoint()
    query_id, query_seq = descs[0], seqs[0]
    duplicates = find_duplicates(descs, seqs)
    descs_merge, seqs_merge = merge_overlaps_and_update(duplicates, query_seq, descs, seqs)
    if genome:
        descs_merge_hs, seqs_merge_hs = highest_score_per_genome_update(descs_merge, seqs_merge, query_seq)
        gap_rem_descs, gap_rem_seqs = gapfilter(descs_merge_hs, seqs_merge_hs, query_seq, gap_ratio_tr)
    else:
        gap_rem_descs, gap_rem_seqs = gapfilter(descs_merge, seqs_merge, query_seq, gap_ratio_tr)
    return gap_rem_descs, gap_rem_seqs

def read_AF_msas(protein, af_output_dir):
    # this inputs the parent directory of alphafold multimer output and finds the location of the MSAs with the 
    # given protein name and returns a list of the msa files. Returns None if didn't find the protein folder
    MSA_dirs = glob.glob(f"{af_output_dir}*/")
    for dir in MSA_dirs:
        with open(f"{dir}msas/chain_id_map.json", 'r') as f:
            data = json.load(f)
        protA, protB = data['A']['description'].split('__')[0], data['B']['description'].split('__')[0]
        if protA == protein:
            return glob.glob(f"{dir}msas/A/*")
        elif protB == protein:
            return glob.glob(f"{dir}msas/B/*")
        elif protA.split('.')[0] == protein: # added to include vaccinia naming UniprotAccession.Genename
            return glob.glob(f"{dir}msas/A/*")
        elif protB.split('.')[0] == protein: # added to include vaccinia naming UniprotAccession.Genename
            return glob.glob(f"{dir}msas/B/*")
    return None

def un_align_sto_or_fasta(seqs):
    # this inputs a list of .sto aligned sequences, seqs, and outputs the raw sqeuences 
    unaligned_seqs = []
    for seq in seqs:
        unaligned_seqs.append(seq.replace('-', ''))
    return unaligned_seqs

def merge_all_overlaps(duplicates, query_seq):
    resolved = {key: None for key in duplicates.keys()}
    while len(duplicates) > 0:
        for key in resolved:
            try:
                if len(duplicates[key]) == 1:
                    resolved[key] = duplicates[key][0][2]
                    resolved[key] = [resolved[key], False] #add tag for not loaded back into MSA
                    del duplicates[key]
            except KeyError:
                continue
        duplicates = resolve_pairwise_overlaps(duplicates, query_seq)
    return resolved

def phmmer_stats_resolved_duplicates(seq_ids, seqs, query_seq, resolved):
#This inputs the list of seq_IDs, the sequences, and the resolved list of duplicate hits from merge_all_overlaps
# It builds a dataframe only adding the resolved sequence for duplicates and calculates the desired stats on the sequences in the df
    query_len = float(len(query_seq.replace('.','').replace('-','')))
    df = pd.DataFrame(columns=['Hits', 'alignments', '%ID', 'coverage', 'gap_ratio'])
    for i in range(0, len(seq_ids)): #iterate through names
        if (seq_ids[i] in resolved) and (resolved[seq_ids[i]][1] == True):
            continue # skip this entry if duplicate has been added to dataframe already
        if (seq_ids[i] in resolved) and (resolved[seq_ids[i]][1] == False):
            seqs[i] = resolved[seq_ids[i]][0]
            resolved[seq_ids[i]][1] = True
        df.loc[i] = [seq_ids[i], seqs[i], ID_score(seqs[i], query_seq), coverage_score(seqs[i], query_len), gap_ratio_score(seqs[i], query_len)]
    return df

def phmmer_sim_filter(df, ID_tr, cov_tr, gap_ratio_tr):
# pulls rows that have values equal to or gretor than %ID and coverage thresholds, and less than gap ratio threshold
    filtered_df = df[(df['%ID']>ID_tr) & (df['coverage']>cov_tr) & (df['gap_ratio']<gap_ratio_tr)]
    rows = len(filtered_df.axes[0])
    return filtered_df, rows


def seed_sto_out(df, file_output):
#this inputs a pd DataFrame with 'Hits' column, 'alisngments' column and saves a STOCKHOLM 1.0 format file
    with open(file_output, 'w') as handle:
        handle.write("# STOCKHOLM 1.0\n")
        max_hit_length = df['Hits'].apply(len).max()
        padding = 4
        for index, row in df.iterrows():
            hit = row['Hits']
            alignment = row['alignments']
            # Calculate the number of spaces needed to align this entry
            num_spaces = max_hit_length - len(hit) + padding
            handle.write(f"{hit}{' ' * num_spaces}{alignment}\n")
        # Write the Stockholm footer
        handle.write("//\n")

def sys_gen_seed_hmmalign(query_fasta, descs, seqs, temp='home/jacob.fenster/temp/', phmmer='phmmer', hmmbuild='hmmbuild', hmmalign='hmmalign', reformat='hhsuite/scripts/remformat.pl'):
    # this inputs a file path to the query fasta, a list of the sequences and descriptions UNALIGNED in fasta
    # conducts a phmmer search of the query against the sequences, generates a seed alignment, converts the seed to a hmm
    # and then aligns the whole list of sequences with the hmm with hmmalign
    # os.system is used for the commands and the paths to the hmmsuite software is input
    # returns list of clean aligned descriptions/sequences in a3m format
    os.makedirs(temp, exist_ok=True) 
    query_descs, query_seqs = bpar.read_fasta(query_fasta)
    query_desc, query_seq = query_descs[0], query_seqs[0]
    # write fasta file of descs and seqs as db for phmmer search
    bpar.fasta_out_single_line(descs, seqs, f"{temp}{query_desc}-unaligned.fasta")
    # generate phmmer search of unaligned sequences
    os.system(f"{phmmer} --notextw -o {temp}{query_desc}-phmmer.log -A {temp}{query_desc}-phmmer.sto {query_fasta} {temp}{query_desc}-unaligned.fasta") 
    # generate the seed alignments from this search 
    sto_descs, sto_seqs = bpar.read_hmmer_sto(f"{temp}{query_desc}-phmmer.sto") 
    phmmer_query_seq = sto_seqs[0]
    duplicates = find_duplicates(sto_descs, sto_seqs)
    resolved = merge_all_overlaps(duplicates, phmmer_query_seq)
    df = phmmer_stats_resolved_duplicates(sto_descs, sto_seqs, phmmer_query_seq, resolved)
    count = len(descs) # get length of query to determine ratio of seed
    # These are threholds for picking a seed sequence determined in Humphreys et al. 2021 Euk Complex prediction
    # pick the strictest set of IDs and gaps to find 25% or more of the original sequences 
    coverage = 0 #they did not filter by this metric
    IDs = [.55, .4, .25]
    gaps = [.2, .35, .5]
    for ID, gap in zip(IDs, gaps):
        seed, rows =  phmmer_sim_filter(df, ID, coverage, gap)
        if rows >= .25*count or rows >= 2500: #take the strictest setting that returns 25% of total or 2500 seq, otherwise go with loosest
            break
    # output the seed alignment 
    seed_sto_out(df, f"{temp}{query_desc}-seed.sto")
    # convert seed alignment to hmm
    os.system(f"{hmmbuild} -o {temp}{query_desc}-hmmbuild.log {temp}{query_desc}-seed.hmm {temp}{query_desc}-seed.sto")
    # align with seed hmm
    os.system(f"{hmmalign} -o {temp}{query_desc}-hmmalign.sto {temp}{query_desc}-seed.hmm {temp}{query_desc}-unaligned.fasta")
    # clean up the alginment 
    aln_descs, aln_seqs = bpar.read_hmmer_sto(f"{temp}{query_desc}-hmmalign.sto")
    aln_query_seq = aln_seqs[0]
    duplicates = find_duplicates(aln_descs, aln_seqs)
    resolved = merge_all_overlaps(duplicates, aln_query_seq)
    clean_aln_descs, clean_aln_seqs = gap_remove_and_filter(aln_descs, aln_seqs, aln_query_seq, resolved, gap_ratio_tr=0.5)
    # convert to a3m format
    bpar.fasta_out_single_line(clean_aln_descs, clean_aln_seqs, f"{temp}{query_desc}-combined.fasta")
    os.system(f"{reformat} fas a3m {temp}{query_desc}-combined.fasta {temp}{query_desc}-combined.a3m")
    final_descs, final_seqs = bpar.read_fasta(f"{temp}{query_desc}-combined.a3m")
    # clean up temp files 
    #os.system(f"rm -r {temp}")
    return final_descs, final_seqs

    
def combine_AF_MSAs(protein_name, af_data_parentdir, temp='home/jacob.fenster/temp', reformat='hhsuite/scripts/reformat.pl'):
    # this inputs a protein name and searches the af_data_parentdir for msas that are identical to that name
    # it then pulls the alignments from each, cleans the data, unalignes the data, and outputs to new lists
    os.makedirs(temp, exist_ok=True)
    master_descs, master_seqs = [], []
    af_msas = read_AF_msas(protein_name, af_data_parentdir)
    if af_msas is None:
        print(f"No AFmult msas for {protein_name}. skipping...")
        return None, None
    else:
        # combine all MSAs
        for msa in af_msas:
            if os.path.basename(msa).split('.')[1] == 'sto':
                if os.path.basename(msa).split('.')[0] == 'pdb_hits': # the pdb hits don't have query sequences included. skipping
                    continue
                temp_descs, temp_seqs = cleanup_sto_or_fasta(msa, gap_ratio_tr=0.5, genome=False)
                temp_seqs = un_align_sto_or_fasta(temp_seqs)
                master_descs.extend(temp_descs), master_seqs.extend(temp_seqs)
            elif os.path.basename(msa).split('.')[1] == 'a3m':
                basename = os.path.basename(msa).split('.')[0]
                new_fasta = f"{temp}{basename}-{protein_name}.fasta"
                os.system(f"{reformat} -v 0 a3m fas {msa} {new_fasta}")
                temp_descs, temp_seqs = cleanup_sto_or_fasta(new_fasta, gap_ratio_tr=0.5, genome=False)
                temp_seqs = un_align_sto_or_fasta(temp_seqs)
                master_descs.extend(temp_descs), master_seqs.extend(temp_seqs)
    #os.system(f"rm -r {temp}")
    master_descs, master_seqs = remove_duplicate_descs(master_descs, master_seqs)   
    return master_descs, master_seqs

def read_a3m_parse_label(fn):
    '''parse an a3m files as a lists of unique labels, non-unique TaxID(codes), and alignments'''
    '''accepts OrthoDB formatted descriptions'''
    '''taken from RoseTTAfold commons originally. Can be edited to exclude certain entries'''
    seq, lab, code = [], [], [] #seq is alignment, lab is unique gene identifier, codes is taxID non-unique identifier
    sequence = '' #placeholder for muli-line fasta 
    is_first = True
    for line in open(fn, "r"):
        if line[0] == '>':
            label = line.strip()[1:]
            #is_incl = True This can be used to exclude certain entries
            if is_first: # include first sequence (query), include full name as code
                is_first = False
                lab.append(label)
                code.append(label)
                continue
            elif "lcl|" in label: #logic to extract genome infomration from NCBI fasta_cds_aa from nuccore
                seq.append(sequence)
                sequence = ''
                lab.append(f"{label.split('|')[1].split('.')[0]}.{label.split('prot_')[1].split('.')[0]}") #>lcl|LT993228.1_prot_SPN68174.1_96/38-763 #endgoal is genome.protein__genome.protein
                code.append(label.split('|')[1].split('.')[0])
            else: #this is the placeholder catch all logic for all other data types. these won't be paired 
                seq.append(sequence) #append previous line and reset
                sequence = ''
                lab.append(label)
                code.append(label) #this is pulling the gene file name
        else:
            sequence += line.rstrip()
    if sequence: #deal with last line 
        seq.append(sequence)
    return seq, lab, code

def dict_builder_longer(seqs, labs, codes):
#this inputs three lists of a3m alignments, unique labels, and non-unique codes and builds a dictionary
#if multiple codes, takes longest of the entries so each code is unique {code:[lab, seq]}
#option to add more logic later for alignments that have multiple hits etc
    result_dict = {}
    for code, lab, seq in zip(codes, labs, seqs):
        if code not in result_dict:
            result_dict[code] = [lab, seq]
        else: #see if current sequence is longer than existing sequence
            existing_seq = result_dict[code][1]
            if len(seq.replace('-','').replace('.','')) > len(existing_seq.replace('-','').replace('.','')): 
                result_dict[code] = [lab, seq]
    return result_dict

def pMSA_dict(dict1, dict2, query1, query2):
#builds a pMSA dict from two dictionaries of single MSAs {code:[lab,seq]}. Paires based on common keys
#concatinates lab for new lab, and seqs for new seq
    pMSA, unique1, unique2 = {}, {}, {}
    # generate query seq!
    pMSA[f"{query1}__{query2}"] = [f"{query1}__{query2}", dict1[query1][1]+dict2[query2][1]]
    for code in dict1:
        if code == query1:
            continue
        elif code in dict2:
            pMSA[code] = [dict1[code][0]+'__'+dict2[code][0], dict1[code][1]+dict2[code][1]]
        else:
            unique1[code] = [dict1[code][0], dict1[code][1]]
    for code in dict2:
        if code == query2: 
            continue
        elif code not in dict1:
            unique2[code] = [dict2[code][0], dict2[code][1]]
    return pMSA, unique1, unique2

def msa_dict_to_a3m(pMSA, output_file):
#this converts a pMSA dictionary {code:[lab,seq]} to an .a3m output with the given basename
#can add option to include unique entries from each MSA
    with open(output_file, 'w') as f:
        for key, value in pMSA.items():
            lab, seq = value
            f.write(f">{lab}\n{seq}\n")
    return
def assemble_paired_unpaired_pMSA(pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs, protA_len, protB_len):
    # this assembles a paired unpaired pMSA based on a3m format input. Adds '_' to the protA and protB seqs based on protA and protB len
    master_descs, master_seqs = [], []
    master_descs.extend(pmsa_descs), master_seqs.extend(pmsa_seqs)
    for desc, seq in zip(protA_descs, protA_seqs):
        master_descs.append(desc)
        master_seqs.append(f"{seq}{'-' * protB_len}")
    for desc, seq in zip(protB_descs, protB_seqs):
        master_descs.append(desc)
        master_seqs.append(f"{'-' * protA_len}{seq}")
    return master_descs, master_seqs

def clean_NCBI_unpaired(pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs, protA_len):
    # this one checks if there are any identical sequences in the unpaired alignments versus the paired alignments 
    
    # de pair
    depairedA, depairedB = [], []
    protA_desc_out, protA_seq_out, protB_desc_out, protB_seq_out = [], [], [], []
    for seq in pmsa_seqs:
        depairedA.append(seq[:protA_len]), depairedB.append(seq[protA_len:])
    for desc, seq in zip(protA_descs, protA_seqs):
        duplicate = False
        for pseq in depairedA:
            if seq == pseq:
                duplicate = True
        if not duplicate:
            protA_desc_out.append(desc), protA_seq_out.append(seq)
    for desc, seq in zip(protB_descs, protB_seqs):
        duplicate = False
        for pseq in depairedB:
            if seq == pseq:
                duplicate = True
        if not duplicate:
            protB_desc_out.append(desc), protB_seq_out.append(seq)
    return protA_desc_out, protA_seq_out, protB_desc_out, protB_seq_out 

def homodimer_MSA(descs, seqs, output_file):
    # this inputs the descs and seqs and writes to the outputfile
    # inserts '-'s to make an unpaired only MSA for homodimer AF runs
    prot_len = len(seqs[0])
    with open(output_file, 'w') as f:
        # generate query homodimer
        f.write(f">{descs[0]}__{descs[0]}\n{seqs[0]}{seqs[0]}\n")
        for desc, seq in zip(descs, seqs):
            f.write(f">{desc}__A\n{seq}{'-' * prot_len}\n")
        for desc, seq in zip(descs, seqs):
            f.write(f">{desc}__B\n{'-' * prot_len}{seq}\n")

def generate_pMSA(proteinA_msa, proteinB_msa, output_dir="", paired_tr=99, unpaired_tr=90,homodimer_tr = 95, temp='/home/jacob.fenster/tmp/',  hhfilter='hhsuite/bin/hhfilter'):
    # given an a3m aligned msa of protein A and B, this generates pMSAs based on NCBI's genome accession
    # and the rest are left unpaired with corresponding number of gaps. 
    # the paired and unpaired sequences are handled separately for hhfilter steps. 
    os.makedirs(temp, exist_ok=True), os.makedirs(f"{temp}pmsa/", exist_ok=True)
    seqA, labA, codeA = read_a3m_parse_label(proteinA_msa)
    dictA = dict_builder_longer(seqA, labA, codeA)
    seqB, labB, codeB = read_a3m_parse_label(proteinB_msa)
    dictB = dict_builder_longer(seqB, labB, codeB)
    protA_len, protB_len = str(len(seqA[0])), str(len(seqB[0]))
    protA_name, protB_name = labA[0], labB[0] # assumes first sequece is query which is the protein name
    pmsa_name = f"{protA_name}-AA{protA_len}__{protB_name}-AA{protB_len}"
    print(f"pMSA generation of {pmsa_name} started")
    # heterodimer logic
    if protA_name != protB_name: 
        pmsa, uniqueA, uniqueB = pMSA_dict(dictA, dictB, protA_name, protB_name)
        # filter each dictionary
        msa_dict_to_a3m(pmsa, f"{temp}pmsa/{pmsa_name}.a3m"),  msa_dict_to_a3m(uniqueA, f"{temp}pmsa/{protA_name}-unique.a3m"), msa_dict_to_a3m(uniqueB, f"{temp}pmsa/{protB_name}-unique.a3m")
        os.system(f"{hhfilter} -v 0 -id {paired_tr} -i {temp}pmsa/{pmsa_name}.a3m -o {temp}pmsa/{pmsa_name}-fil.a3m")
        os.system(f"{hhfilter} -v 0 -id {unpaired_tr} -i {temp}pmsa/{protA_name}-unique.a3m -o {temp}pmsa/{protA_name}-fil.a3m")
        os.system(f"{hhfilter} -v 0 -id {unpaired_tr} -i {temp}pmsa/{protB_name}-unique.a3m -o {temp}pmsa/{protB_name}-fil.a3m")
        # assemble paired-unpaired pMSA. Set to empty lists if no file (ie no seqs) 
        try:
            pmsa_descs, pmsa_seqs = bpar.read_fasta(f"{temp}pmsa/{pmsa_name}-fil.a3m")
        except FileNotFoundError:
            pmsa_descs, pmsa_seqs = [], []
        try:
            protA_descs, protA_seqs = bpar.read_fasta(f"{temp}pmsa/{protA_name}-fil.a3m")
        except FileNotFoundError:
            protA_descs, protA_seqs = [], []
        try:
            protB_descs, protB_seqs = bpar.read_fasta(f"{temp}pmsa/{protB_name}-fil.a3m")
        except FileNotFoundError:
            protB_descs, protB_seqs = [], []
        protA_descs_clean, protA_seqs_clean, protB_descs_clean, protB_seqs_clean = clean_NCBI_unpaired(pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs, int(protA_len)) 
        pu_descs, pu_seqs = assemble_paired_unpaired_pMSA(pmsa_descs, pmsa_seqs, protA_descs_clean, protA_seqs_clean, protB_descs_clean, protB_seqs_clean, int(protA_len), int(protB_len))
        # ouptut paired and paired_unpaired pMSAs 
        os.makedirs(f"{output_dir}paired/", exist_ok=True), os.makedirs(f"{output_dir}paired_unpaired/", exist_ok=True)
        bpar.fasta_out_single_line(pu_descs, pu_seqs, f"{output_dir}paired_unpaired/{pmsa_name}-pu.a3m")
        bpar.fasta_out_single_line(pmsa_descs, pmsa_seqs, f"{output_dir}paired/{pmsa_name}-p.a3m")
        #return pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs
    # homodimer logic, only generate unpaired msa
    elif protA_name == protB_name: 
        os.system(f"{hhfilter} -v 0 -id {homodimer_tr} -i {proteinA_msa} -o {temp}pmsa/{protA_name}-fil.a3m")
        try:
            protA_descs, protA_seqs = bpar.read_fasta(f"{temp}pmsa/{protA_name}-fil.a3m")
        except FileNotFoundError:
            protA_descs, protA_seqs = [], []
        os.makedirs(f"{output_dir}homodimer/", exist_ok=True)
        homodimer_MSA(protA_descs, protA_seqs, f"{output_dir}homodimer/{pmsa_name}-h.a3m")
    print(f"pMSA generation of {pmsa_name} complete")
        #return 

def sub_generate_pMSA(proteinA_msa, proteinB_msa, output_dir="", paired_tr=99, unpaired_tr=90,homodimer_tr = 95, temp='/home/jacob.fenster/tmp/',  hhfilter='hhsuite/bin/hhfilter'):
    # this one uses subprcoess to record stdout and error text
    # given an a3m aligned msa of protein A and B, this generates pMSAs based on NCBI's genome accession
    # and the rest are left unpaired with corresponding number of gaps. 
    # the paired and unpaired sequences are handled separately for hhfilter steps. 
    os.makedirs(temp, exist_ok=True), os.makedirs(f"{temp}pmsa/", exist_ok=True)
    seqA, labA, codeA = read_a3m_parse_label(proteinA_msa)
    dictA = dict_builder_longer(seqA, labA, codeA)
    seqB, labB, codeB = read_a3m_parse_label(proteinB_msa)
    dictB = dict_builder_longer(seqB, labB, codeB)
    protA_len, protB_len = str(len(seqA[0])), str(len(seqB[0]))
    protA_name, protB_name = labA[0], labB[0] # assumes first sequece is query which is the protein name
    if '.' in protA_name or '.' in protB_name:
        pmsa_name = f"{protA_name.split('.')[0]}-AA{protA_len}__{protB_name.split('.')[0]}-AA{protB_len}"
    else:
        pmsa_name = f"{protA_name}-AA{protA_len}__{protB_name}-AA{protB_len}"
    print(f"pMSA generation of {pmsa_name} started")
    # heterodimer logic
    if protA_name != protB_name: 
        pmsa, uniqueA, uniqueB = pMSA_dict(dictA, dictB, protA_name, protB_name)
        # filter each dictionary
        msa_dict_to_a3m(pmsa, f"{temp}pmsa/{pmsa_name}.a3m"),  msa_dict_to_a3m(uniqueA, f"{temp}pmsa/{protA_name}-unique.a3m"), msa_dict_to_a3m(uniqueB, f"{temp}pmsa/{protB_name}-unique.a3m")
        pmsa_filter = subprocess.run(f"{hhfilter} -v 0 -id {paired_tr} -i {temp}pmsa/{pmsa_name}.a3m -o {temp}pmsa/{pmsa_name}-fil.a3m", shell=True, text=True, capture_output=True)
        protA_filter = subprocess.run(f"{hhfilter} -v 0 -id {unpaired_tr} -i {temp}pmsa/{protA_name}-unique.a3m -o {temp}pmsa/{protA_name}-fil.a3m", shell=True, capture_output=True)
        protB_filter = subprocess.run(f"{hhfilter} -v 0 -id {unpaired_tr} -i {temp}pmsa/{protB_name}-unique.a3m -o {temp}pmsa/{protB_name}-fil.a3m", shell=True, capture_output=True)
        print(f"{pmsa_name} stdout: {pmsa_filter.stdout}; stderr: {pmsa_filter.stderr}")
        print(f"{protA_name} stdout: {protA_filter.stdout}; stderr: {protA_filter.stderr}")
        print(f"{protB_name} stdout: {protB_filter.stdout}; stderr: {protB_filter.stderr}")
        # assemble paired-unpaired pMSA. Set to empty lists if no file (ie no seqs) 
        try:
            pmsa_descs, pmsa_seqs = bpar.read_fasta(f"{temp}pmsa/{pmsa_name}-fil.a3m")
        except FileNotFoundError:
            pmsa_descs, pmsa_seqs = [], []
        try:
            protA_descs, protA_seqs = bpar.read_fasta(f"{temp}pmsa/{protA_name}-fil.a3m")
        except FileNotFoundError:
            protA_descs, protA_seqs = [], []
        try:
            protB_descs, protB_seqs = bpar.read_fasta(f"{temp}pmsa/{protB_name}-fil.a3m")
        except FileNotFoundError:
            protB_descs, protB_seqs = [], []
        protA_descs_clean, protA_seqs_clean, protB_descs_clean, protB_seqs_clean = clean_NCBI_unpaired(pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs, int(protA_len)) 
        pu_descs, pu_seqs = assemble_paired_unpaired_pMSA(pmsa_descs, pmsa_seqs, protA_descs_clean, protA_seqs_clean, protB_descs_clean, protB_seqs_clean, int(protA_len), int(protB_len))
        # ouptut paired and paired_unpaired pMSAs 
        os.makedirs(f"{output_dir}paired/", exist_ok=True), os.makedirs(f"{output_dir}paired_unpaired/", exist_ok=True)
        bpar.fasta_out_single_line(pu_descs, pu_seqs, f"{output_dir}paired_unpaired/{pmsa_name}-pu.a3m")
        bpar.fasta_out_single_line(pmsa_descs, pmsa_seqs, f"{output_dir}paired/{pmsa_name}-p.a3m")
        #return pmsa_descs, pmsa_seqs, protA_descs, protA_seqs, protB_descs, protB_seqs
    # homodimer logic, only generate unpaired msa
    elif protA_name == protB_name: 
        protA_filter = subprocess.run(f"{hhfilter} -v 0 -id {homodimer_tr} -i {proteinA_msa} -o {temp}pmsa/{protA_name}-fil.a3m", shell=True, capture_output=True)
        print(f"{protA_name} stdout: {protA_filter.stdout}; stderr: {protA_filter.stderr}")
        try:
            protA_descs, protA_seqs = bpar.read_fasta(f"{temp}pmsa/{protA_name}-fil.a3m")
        except FileNotFoundError:
            protA_descs, protA_seqs = [], []
        os.makedirs(f"{output_dir}homodimer/", exist_ok=True)
        homodimer_MSA(protA_descs, protA_seqs, f"{output_dir}homodimer/{pmsa_name}-h.a3m")
    print(f"pMSA generation of {pmsa_name} complete")

def combined_msa(protein_name, combined_msa_out, jhmmer_dir, query_fasta_dir, afmult_out, temp, reformat):
    # generates a combined MSA. Def to be parallelized 
    print(f"starting combined MSA for {protein_name}")
    jhmmer_file = f"{jhmmer_dir}{protein_name}_jhmmer.sto"
    query_fasta = f"{query_fasta_dir}{protein_name}.fasta"
    master_seqs, master_descs = [], []
    viral_descs, viral_seqs = cleanup_sto_or_fasta(jhmmer_file, gap_ratio_tr=0.5, genome=True)
    print(f"read jackhmmer sto file")
    if viral_descs[0] == 'query': # fix old jhmmer files that had 'query' as query
        viral_descs[0] = protein_name
    viral_seqs = un_align_sto_or_fasta(viral_seqs) 
    master_descs.extend(viral_descs), master_seqs.extend(viral_seqs)
    # read and combine  MSAs from AFmult output
    af_descs, af_seqs = combine_AF_MSAs(protein_name, afmult_out, temp=temp, reformat=reformat)
    print("read AF MSAs")
    if af_descs is None:
        print(f"No AFmult msas for {protein_name}. skipping...")
    else:
        master_descs.extend(af_descs), master_seqs.extend(af_seqs)
        aln_descs, aln_seqs = sys_gen_seed_hmmalign(query_fasta, master_descs, master_seqs, temp=f"{temp}aln/", phmmer='phmmer', hmmbuild='hmmbuild', hmmalign='hmmalign', reformat=reformat)
        print("generated alignment")
        # remove duplicate alignments prioritizing those in the all NCBI virus db
        unique_aln_descs, unique_aln_seqs = remove_duplicates_prioitize_NCBI(aln_descs, aln_seqs)
        # send a3m MSA to output file  
        bpar.fasta_out_single_line(unique_aln_descs, unique_aln_seqs, f"{combined_msa_out}{protein_name}.a3m")
        print(f"Combined MSA for {protein_name} complete")
