import sys
import os
from collections import defaultdict

def read_hmmer_sto(input_sto):
#this reads a Stockholm 1.0 format MSA and parses between blocks and allows for redundant hits
    query = True #holder to find query
    query_id, query_seq = '', ''
    seq_ids = [] 
    #first extract the query seqeunce and id list of all alignments
    with open(input_sto, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('#') or line.startswith('\n'):
            if line.startswith('#=GC RF'):
                break #stop after first block, should work if only one block
            continue
        else:
            seq_id, seq = line.strip().split(None, 1)
            seq_ids.append(seq_id)
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

def find_alifrom_alito(seq):
    alnlen = len(seq)
    alifrom, alito = None, None #index of start of alignment and end of alignment. First index of alignment is zero
    for c in range(alnlen): #iterate through alignment columns
        if seq[c].isalpha(): #find first entry of residue
            if alifrom is None:
                alifrom = c
            alito = c #pulls the last alpha character in the list
    return alifrom, alito

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

def overlapping_ID_score(seq, overlap, query_seq):
    #scores the percent ID of the given overlapping region. Input is full alignment of seq and query. this slices overlap
    seq, query_seq = seq[overlap[0]:overlap[1]+1], query_seq[overlap[0]:overlap[1]+1]
    identical, aligned = 0, 0
    for c in range(len(query_seq)):
        if seq[c].isupper():
            aligned +=1
            if seq[c] == query_seq[c]:
                identical += 1
    ID = identical/aligned
    return ID

def gap_ratio_score(seq, query_seq):
#scores the gap ratio of seq versus query seq. input is full length sequence and scores full lenth
#input is query GAP column REMOVED sequences
    aln_len = len(query_seq)
    gap = 0
    for c in range(aln_len):
        if seq[c] == '-':
            gap += 1
    gap_ratio = gap/aln_len
    return gap_ratio


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


def gap_remove_and_filter_FASTAout(seq_ids, seqs, query_seq, resolved, output_filename, gap_ratio_tr):
    gap_indicies = [i for i, char in enumerate(query_seq) if char in ['-']] #remove all - columns 
    gap_rem_query = ''.join([char for i, char in enumerate(query_seq) if i not in gap_indicies])
    with open(output_filename, 'w') as handle:
        for seq_id, seq in zip(seq_ids, seqs):
            if seq_id in resolved:
                if not resolved[seq_id][1]:
                    gap_rem_seq = ''.join([char for i, char in enumerate(resolved[seq_id][0]) if i not in gap_indicies])
                    if gap_ratio_score(gap_rem_seq, gap_rem_query) < gap_ratio_tr:
                        handle.write(f">{seq_id}\n{gap_rem_seq}\n")
                    resolved[seq_id][1] = True 
                continue
            gap_rem_seq = ''.join([char for i, char in enumerate(seq) if i not in gap_indicies])
            if gap_ratio_score(gap_rem_seq, gap_rem_query) < gap_ratio_tr:
                handle.write(f">{seq_id}\n{gap_rem_seq}\n")
    return

#usage: $ python3 hmmalign_cleanup.py <input.sto alignment> <gap ratio threshold, 0.5> <resolved fasta format filename>

input_sto, gap_ratio_tr, output_filename = sys.argv[1], float(sys.argv[2]), sys.argv[3]
seq_ids, seqs = read_hmmer_sto(input_sto)
query_seq = seqs[0] #assumes the first alignment is the query sequence
duplicates = find_duplicates(seq_ids, seqs)
resolved = merge_all_overlaps(duplicates, query_seq)
gap_remove_and_filter_FASTAout(seq_ids, seqs, query_seq, resolved, output_filename, gap_ratio_tr)