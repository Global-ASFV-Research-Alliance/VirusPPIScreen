import os
import pandas as pd

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

def parse_uniprot_proteome(proteome_file):
    # this parses the accession number for each entry in a UniProt proteome file and returns a list
    proteome_df = pd.DataFrame(columns=['description', 'sequence'])
    proteome_df.index.name = 'uniprot_accession'
    descs, seqs = read_fasta(proteome_file)
    for desc, seq in zip(descs, seqs):
       accession = desc.split('|')[1]
       proteome_df.loc[accession] = {'description':desc, 'sequence':seq}
    return proteome_df 

