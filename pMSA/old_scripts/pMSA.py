import os
import sys
import itertools
import numpy as np

def get_query_len(fn):
#calculates the length of the query, first entry in .a3m or fasta
    stop = False
    seq, seqlen = '', ''
    for line in open(fn, "r"): #calculate length of query for both files
        if line[0] == '>':
            if stop == True:
                seqlen = len(seq)
                break
            stop = True
            continue
        seq+=(line.rstrip())
    return seqlen

def read_a3m(fn):
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
            elif "pub_og_id" in label: #logic for if OrthoDB MSA
                seq.append(sequence) #append previous line and reset
                sequence = ''
                lab.append(label.split(' ')[0]) #use full TaxID:geneID as unique label
                code.append(label.split(':')[0]) #use OrthoDB taxID as non-unique code. fmt is taxID_ver#
            elif ":" in label: #logic for OrthoDB -> hmmalign-ed descriptions ex. ">1128424_1:000067\n"
                seq.append(sequence) #append previous line and reset
                sequence = ''
                lab.append(label)
                code.append(label.split(':')[0])
            else: #this is the placeholder catch all logic for FastaBank internal ASFV data
                seq.append(sequence) #append previous line and reset
                sequence = ''
                lab.append(label+'_'+os.path.basename(fn).split('.')[0].split('_')[0])
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

def pMSA_dict(dict1, dict2):
#builds a pMSA dict from two dictionaries of single MSAs {code:[lab,seq]}. Paires based on common keys
#concatinates lab for new lab, and seqs for new seq
    pMSA, unique1, unique2 = {}, {}, {}
    for code in dict1:
        if code in dict2:
            pMSA[code] = [dict1[code][0]+'__'+dict2[code][0], dict1[code][1]+dict2[code][1]]
        else:
            unique1[code] = [dict1[code][0], dict1[code][1]]
    for code in dict2:
        if code not in dict1:
            unique2[code] = [dict2[code][0], dict2[code][1]]
    return pMSA, unique1, unique2

def pMSA_dict_to_a3m(pMSA, paired_basename):
#this converts a pMSA dictionary {code:[lab,seq]} to an .a3m output with the given basename
#can add option to include unique entries from each MSA
    with open(paired_basename+'-pMSA.a3m', 'w') as f:
        for key, value in pMSA.items():
            lab, seq = value
            f.write(f">{lab}\n{seq}\n")
    return

MSA_files = sys.argv[1:] #command line list of files. Can be changed to accept a .csv
pairs = list(itertools.combinations(MSA_files, 2)) #all pairwise combinations of input files
for pair in pairs:
    querylen1, querylen2 = str(get_query_len(pair[0])), str(get_query_len(pair[1]))
    seq1, lab1, code1 = read_a3m(pair[0])
    dict1 = dict_builder_longer(seq1, lab1, code1)
    seq2, lab2, code2 = read_a3m(pair[1])
    dict2 = dict_builder_longer(seq2, lab2, code2)
    pMSA, unique1, unique2 = pMSA_dict(dict1, dict2)
    #determine output basename
    if '10239' in pair[0]: #check to see if orthoDB ASFV level. specific to my nonsense naming scheme
        name1 = os.path.basename(pair[0]).split('__')[0].split('_')[1]
        name2 = os.path.basename(pair[1]).split('__')[0].split('_')[1]
        paired_basename = name1+'-AA'+querylen1+'__'+name2+'-AA'+querylen2+'-'+ os.path.basename(pair[0]).split('.')[0].split('10239')[1].split('_')[0]
    elif '_rearranged_clustalo_filtered' in pair[0]:
        name1=os.path.basename(pair[0]).split('.')[0].split('_rearranged_clustalo_filtered')[0]
        name2=os.path.basename(pair[1]).split('.')[0].split('_rearranged_clustalo_filtered')[0]
        paired_basename = name1+'-AA'+querylen1+'__'+name2+'-AA'+querylen2
    else: #takes full basename as gene/query name
        name1 = os.path.basename(pair[0]).split('.')[0]
        name2 = os.path.basename(pair[1]).split('.')[0]
        paired_basename = f"{name1}-AA{querylen1}__{name2}-AA{querylen2}"
    pMSA_dict_to_a3m(pMSA, paired_basename)





