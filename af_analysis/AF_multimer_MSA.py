import os, sys, glob, json, shutil
import pdb
# This script moves precomputed MSAs into new AFmultimer output directories
# useful for when you have the monomers computed and want to compute new combinations of PPIs
# right now only works for heterodimers, not homodimers

#this script is written for the given format of input fasta multimer peptides
# input fasta is geneA__geneB.fasta with descs:
# >geneA
# >geneB
# and the input MSAs are in directory geneA__geneB/
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

in_MSA_dir, output_dir, fasta = sys.argv[1], sys.argv[2], sys.argv[3]
descs, seqs = read_fasta(fasta)
MSA_dirs = glob.glob(f"{in_MSA_dir}*/")
fasta_basename = os.path.basename(fasta).split('.')[0]
# make output directories
try:
    os.makedirs(f"{output_dir}{fasta_basename}/")
    os.makedirs(f"{output_dir}{fasta_basename}/msas/")
except FileExistsError:
    # see if the msas have been moved successfully but af was not run
    if os.path.exists(f"{output_dir}{fasta_basename}/ranked_0.pdb"):
        print(f"output solution already exists for {fasta_basename}")
        sys.exit(1)
    elif os.path.exists(f"{output_dir}{fasta_basename}/msas/A") and os.path.exists(f"{output_dir}{fasta_basename}/msas/B"):
        print(f"output MSAs already exists for {fasta_basename}, continue with calculation")
        sys.exit(0)
    else:
        print(f"some other error occurred for {fasta_basename}, skipping")
        sys.exit(1)

for dir in MSA_dirs:
    with open(f"{dir}msas/chain_id_map.json", 'r') as f:
        data = json.load(f)
    protA, protB = data['A']['description'], data['B']['description']
    if protA.split('__')[0] == protB.split('__')[0]: # handle homodimer case
        protA, protB = protA.split('__')[0], protB.split('__')[0]
    if protA == descs[0] and not os.path.exists(f"{output_dir}{fasta_basename}/msas/A/") and os.path.exists(f"{dir}msas/A/"):
        shutil.copytree(f"{dir}msas/A/", f"{output_dir}{fasta_basename}/msas/A/")
    elif protA == descs[1] and not os.path.exists(f"{output_dir}{fasta_basename}/msas/B/") and os.path.exists(f"{dir}msas/A/"):
        shutil.copytree(f"{dir}msas/A/", f"{output_dir}{fasta_basename}/msas/B/")
    if protB == descs[0] and not os.path.exists(f"{output_dir}{fasta_basename}/msas/A/") and os.path.exists(f"{dir}msas/B/"):
        shutil.copytree(f"{dir}msas/B/", f"{output_dir}{fasta_basename}/msas/A/")
    elif protB == descs[1] and not os.path.exists(f"{output_dir}{fasta_basename}/msas/B/") and os.path.exists(f"{dir}msas/B/"):
        shutil.copytree(f"{dir}msas/B/", f"{output_dir}{fasta_basename}/msas/B/")

if not os.path.exists(f"{output_dir}{fasta_basename}/msas/A/"):
    print(f"{descs[0]} msas not found")
    exit_status = 1
if not os.path.exists(f"{output_dir}{fasta_basename}/msas/B/"):
    print(f"{descs[1]} msas not found")
    exit_status = 1
if os.path.exists(f"{output_dir}{fasta_basename}/msas/A/") and os.path.exists(f"{output_dir}{fasta_basename}/msas/B/"):
    print(f"successfully transferred MSAs for {fasta_basename}.fasta")
    exit_status = 0
sys.exit(exit_status)
