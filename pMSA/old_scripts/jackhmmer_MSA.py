import os
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import pdb 

def read_hmmer_sto(input_sto):
#this reads a Stockholm 1.0 format MSA and parses between blocks and allows for redundant hits
# Uses '#=GC RF' as end of block identifier. This is likely HMMER specific 
# jackhmmer outputs the query alignemnt as the first line in the alignments
    query, query_id = True, '' #holder to find query for block parsing
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

def ID_score(seq, query_seq):
    aligned, identical = 0, 0
    for c in range(len(query_seq)):
        if seq[c].isupper():
            aligned += 1
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

def find_duplicates(seq_ids, seqs):
    #record indicies in original list, adding multiple idicies for duplicate hits (same protein)
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
def read_hmmer_sto(input_sto):
#this reads a Stockholm 1.0 format MSA and parses between blocks and allows for redundant hits
# Uses '#=GC RF' as end of block identifier. This is likely HMMER specific 
# jackhmmer outputs the query alignemnt as the first line in the alignments
    query, query_id = True, '' #holder to find query for block parsing
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

def gapfilter_FASTAout(seq_ids, seqs, query_seq, output_filename, gap_ratio_tr):
    gap_indicies = [i for i, char in enumerate(query_seq) if char in ['-']] #remove all - columns 
    gap_rem_query = ''.join([char for i, char in enumerate(query_seq) if i not in gap_indicies])
    with open(output_filename, 'w') as handle:
        for seq_id, seq in zip(seq_ids, seqs):
            gap_rem_seq = ''.join([char for i, char in enumerate(seq) if i not in gap_indicies])
            if gap_ratio_score(gap_rem_seq, gap_rem_query) < gap_ratio_tr:
                handle.write(f">{seq_id}\n{gap_rem_seq}\n")
    return
        



input_sto, gap_ratio_tr, output_filename = sys.argv[1], float(sys.argv[2]), sys.argv[3]
#read the input .sto file to lists 
seq_ids, seqs = read_hmmer_sto(input_sto)
#jackhmmer outputs the first line as the query alignment
query_id, query_seq = seq_ids[0], seqs[0]
#find alignments that target identical proteins
duplicates = find_duplicates(seq_ids, seqs)

#concatenate and merge alignments that target identical proteins. Update original list with merged alignments
seq_ids_merge, seqs_merge = merge_overlaps_and_update(duplicates, query_seq, seq_ids, seqs)

#find instances where there are multiple proteins hits within the same genome but of different proteins
#resolve by higher ID%, removing the other entry
seq_ids_unique, seqs_unique = highest_score_per_genome_update(seq_ids_merge, seqs_merge, query_seq)
gapfilter_FASTAout(seq_ids_unique, seqs_unique, query_seq, output_filename, gap_ratio_tr)
#need to test the algorithm to this point to see if it is functioning as expected 

