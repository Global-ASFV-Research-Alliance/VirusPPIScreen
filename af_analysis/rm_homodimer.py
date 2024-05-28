import os, sys, glob
tranche1 = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/afmult_ian-tophits/leftovers/tranche1/"
tranche2 = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/afmult_ian-tophits/leftovers/tranche2/"

fasta_t1 = glob.glob(f"{tranche1}*.fasta")
fasta_t2 = glob.glob(f"{tranche2}*.fasta")
for fasta in fasta_t1:
    basename = os.path.basename(fasta).split('.')[0]
    protA, protB = basename.split('__')[0], basename.split('__')[1]
    if protA == protB:
        os.system(f"rm {fasta}")
for fasta in fasta_t2:
    basename = os.path.basename(fasta).split('.')[0]
    protA, protB = basename.split('__')[0], basename.split('__')[1]
    if protA == protB:
        os.system(f"rm {fasta}")
