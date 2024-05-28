import os, sys, glob, pdb, time
from collections import defaultdict
import numpy as np
import pandas as pd
import itertools

"""
This script is a combination of jackhmmer_MSA and pMSA where the input data is just jackhmmer searches of the all virus NCBI set
the pMSA algorithm now pairs by genome accession and then keeps all unpaired alignments adding gaps on the unpaired segment
These paired+unpaired pMSAs are then clustered with hhfilter at 99%
Created 20231212
"""

"""
jackhmmer_MSA defs
"""
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

def gapfilter(seq_ids, seqs, query_seq, gap_ratio_tr):
    gapfilter_seqids, gapfilter_seqs = [], []
    gap_indicies = [i for i, char in enumerate(query_seq) if char in ['-']] #remove all - columns 
    gap_rem_query = ''.join([char for i, char in enumerate(query_seq) if i not in gap_indicies])
    for seq_id, seq in zip(seq_ids, seqs):
        gap_rem_seq = ''.join([char for i, char in enumerate(seq) if i not in gap_indicies])
        if gap_ratio_score(gap_rem_seq, gap_rem_query) < gap_ratio_tr:
            gapfilter_seqids.append(seq_id)
            gapfilter_seqs.append(gap_rem_seq)
    return gapfilter_seqids, gapfilter_seqs

"""
pMSA defs
"""

def parse_seq_lab_code(seq_ids, seqs):
    """
    This def inputs the raw seq_ids and seq of a filtered MSA and outputs the parsed seqs label and code for pairing
    right now it is only programmed to accept all viral genomes and parse labels and codes for these
    """
    seq, lab, code = [], [], []
    for i in range(len(seq_ids)):
        if i == 0: # include first sequence (query), include full name as code
            seq.append(seqs[i]), lab.append(seq_ids[i]), code.append(seq_ids[i])
        elif "lcl|" in seq_ids[i]:
            seq.append(seqs[i])
            lab.append(f"{seq_ids[i].split('|')[1].split('.')[0]}.{seq_ids[i].split('prot_')[1].split('.')[0]}") #>lcl|LT993228.1_prot_SPN68174.1_96/38-763 #endgoal is genome.protein__genome.protein
            code.append(seq_ids[i].split('|')[1].split('.')[0])
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

def pMSA_paired_and_unpaired(dict1, dict2, querylen1, querylen2):
#builds a pMSA dict from two dictionaries of single MSAs {code:[lab,seq]}. Paires based on common keys
#Alignments must be in FASTA format
#concatinates lab for new lab, and seqs for new seq
#if label unique to dict1 or dict2 it adds it with the corresponding gaps according to querylen1 and querylen2 length
    pMSA = {}
    for code in dict1:
        if code in dict2: #pair and concatinate
            pMSA[code] = [dict1[code][0]+'__'+dict2[code][0], dict1[code][1]+dict2[code][1]]
        else: # add gaps for peptide 2
            pMSA[code] = [f"{dict1[code][0]}__", f"{dict1[code][1]}{'-' * querylen2}"]
    for code in dict2:
        if code not in dict1: #add gaps for peptide 1
            pMSA[code] = [f"__{dict2[code][0]}", f"{'-' * querylen1}{dict2[code][1]}"]
    return pMSA

def pMSA_dict_to_a3m(pMSA, paired_basename):
#this converts a pMSA dictionary {code:[lab,seq]} to an .a3m output with the given basename
#can add option to include unique entries from each MSA
    with open(paired_basename+'-pMSA.a3m', 'w') as f:
        for key, value in pMSA.items():
            lab, seq = value
            f.write(f">{lab}\n{seq}\n")
    return

def pMSA_dict_to_fasta(pMSA, output_path):
#this converts a pMSA dictionary {code:[lab,seq]} to an .a3m output with the given basename
#can add option to include unique entries from each MSA
    with open(output_path, 'w') as f:
        for key, value in pMSA.items():
            lab, seq = value
            f.write(f">{lab}\n{seq}\n")
    return


"""
Main script
"""
gap_ratio_tr = 0.5
id_clust = 99
input_dir = "/mnt/c/Users/jacob.fenster/USDA/Spinard, Edward - REE-ARS - PPI/JacobsFolder/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/jhmmer_sto"
temp_dir = "/mnt/c/tmp" # to hold raw fasta pMSAs
output_dir = "/mnt/c/Users/jacob.fenster/USDA/Spinard, Edward - REE-ARS - PPI/JacobsFolder/pMSA_algorithm/output/jackhmmer_complete_viral_genomes_NCBI/ASFV_Georgia2007/paired&unpaired_clust99/pMSAs"
reformat = "/mnt/c/Users/jacob.fenster/hhsuite/scripts/reformat.pl"
hhfilter = "/mnt/c/Users/jacob.fenster/hhsuite/bin/hhfilter"

input_stos = glob.glob(f"{input_dir}/*.sto")
MSAs = {} # dict to hold {basename:([ids], [seqs])}
#read each jackhmmer sto, filter, and load onto dictionary
print(f"Paired and Unpaired pMSA started with {id_clust}% clustering\nThere are {len(input_stos)} jackhmmer sto files\ninput sto: {input_dir}\ntemp dir: {temp_dir}\noutput dir: {output_dir}")
timei= time.time()
for sto in input_stos:
    # split the basename to remove the _jhmmer tag
    # specific to this current run
    basename = os.path.basename(sto).split('_jhmmer')[0]
    #read the input .sto file to lists 
    seq_ids, seqs = read_hmmer_sto(sto)
    #jackhmmer outputs the first line as the query alignment
    query_id, query_seq = seq_ids[0], seqs[0]
    #find alignments that target identical proteins
    duplicates = find_duplicates(seq_ids, seqs)
    #concatenate and merge alignments that target identical proteins. Update original list with merged alignments
    seq_ids_merge, seqs_merge = merge_overlaps_and_update(duplicates, query_seq, seq_ids, seqs)
    #find instances where there are multiple proteins hits within the same genome but of different proteins
    #resolve by higher ID%, removing the other entry
    seq_ids_unique, seqs_unique = highest_score_per_genome_update(seq_ids_merge, seqs_merge, query_seq)
    gapfilter_seqids, gapfilter_seqs = gapfilter(seq_ids_unique, seqs_unique, query_seq, gap_ratio_tr)
    MSAs[basename] = (gapfilter_seqids, gapfilter_seqs)
print(f"The input sto filtering is complete")
# build pMSAs. have not made progress past this point. Need to re check 
# MSAs[key][0] is list of seq_ids, MSAs[key][1] is list of alignments
# MSAs[key][0][0] is query ID, MSAs[key][1][0] is query alignment
# the MSA is in fasta format
timefil = time.time()
pairs = list(itertools.combinations(MSAs.keys(), 2)) #all pairwise combinations of input files
for pair in pairs:
    querylen1, querylen2 = len(MSAs[pair[0]][1][0]), len(MSAs[pair[1]][1][0])
    seq1, lab1, code1 = parse_seq_lab_code(MSAs[pair[0]][0], MSAs[pair[0]][1])
    dict1 = dict_builder_longer(seq1, lab1, code1)
    seq2, lab2, code2 = parse_seq_lab_code(MSAs[pair[1]][0], MSAs[pair[1]][1])
    dict2 = dict_builder_longer(seq2, lab2, code2)
    #again, these should be in fasta format for now. need to figure out what to do 
    # also need to add the addition of gaps to the unpaired pMSAs
    pMSA = pMSA_paired_and_unpaired(dict1, dict2, querylen1, querylen2)
    paired_basename = f"{pair[0]}-AA{querylen1}__{pair[1]}-AA{querylen2}"
    output_path = f"{temp_dir}/{paired_basename}.fasta"
    pMSA_dict_to_fasta(pMSA, output_path)
    # convert to a3m format
    os.system(f"{reformat} fas a3m {output_path} {temp_dir}/{paired_basename}.a3m")
    # cluster at 99% ID where pu99 is paired and unpaired 99% ID clustering
    os.system(f"{hhfilter} -v 0 -id {id_clust} -i {temp_dir}/{paired_basename}.a3m -o \"{output_dir}/{paired_basename}-pMSApu99.a3m\"")
timefinal = time.time()
print(f"script complete. Runtime was {(timefil-timei)/60:.1} minutes for filtering and {(timefinal-timefil)/60:.1} minutes for pMSA generation, reformatting and clustering")

