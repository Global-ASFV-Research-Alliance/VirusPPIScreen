
import os
import sys
import pdb

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a list of descs and a list of seqs.
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

def unique_strings_with_descs(alns, descs):
    unique_dict = {}
    for string, description in zip(alns, descs):
        if description == 'query':
            unique_dict[string] = description
        elif string not in unique_dict:
            unique_dict[string] = description
    unique_alns = list(unique_dict.keys())
    unique_descs = list(unique_dict.values())
    return unique_alns, unique_descs

a3m_file, output_file = sys.argv[1], sys.argv[2]
descs, alns = read_fasta(a3m_file)
unique_alns, unique_descs = unique_strings_with_descs(alns, descs)
with open(output_file, 'w') as f:
    for desc, aln in zip(unique_descs, unique_alns):
        f.write(f">{desc}\n{aln}\n")

