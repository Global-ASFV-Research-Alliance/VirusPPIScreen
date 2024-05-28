import os, glob, pdb, itertools
import common.bio_parsers as bpar
# This generates the fasta files for AF homodimer input for MSA generation of FMDV O1 Campos 
fmdv_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/FMDV/FMDV-O1Campos_proteome.fa"
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/alphafold_multimer/homodimers_O1Campos/"
descs, seqs = bpar.read_fasta(fmdv_proteome)
if False:
    for desc, seq in zip(descs, seqs):
        with open(f"{output_dir}{desc}__{desc}.fasta", 'w') as f:
            f.write(f">{desc}__A\n{seq}\n>{desc}__B\n{seq}")
# generate AFmultimer heterodimerset 
protein_names = []
protein_dict = {}
descs, seqs = bpar.read_fasta(fmdv_proteome)
for desc, seq in zip(descs, seqs):
    protein, basename = desc, desc
    protein_names.append(protein)
    protein_dict[desc] = seq

homo_heterodimer_output = "/lustrefs/fadru/projects/asfv-ppi/Input/FMDV/alphafold_multimer/hetero_homodimers-O1Campos/"
pairs = list(itertools.combinations(protein_names, 2)) #all pairwise combinations of input files
for pair in pairs:
    protA, protB = min([pair[0], pair[1]]), max([pair[0], pair[1]])
    ppi = f"{protA}__{protB}"
    with open(f"{homo_heterodimer_output}{ppi}.fasta", 'w') as f:
        f.write(f">{protA}\n{protein_dict[protA]}\n>{protB}\n{protein_dict[protB]}")
# generate AFmultimer homodimerset 
for protein in protein_names:
    ppi = f"{protein}__{protein}"
    with open(f"{homo_heterodimer_output}{ppi}.fasta", 'w') as f:
        f.write(f">{protein}__A\n{protein_dict[protein]}\n>{protein}__B\n{protein_dict[protein]}")
