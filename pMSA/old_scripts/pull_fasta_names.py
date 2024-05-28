import os
import common.bio_parsers as bpar

descs, seqs = bpar.read_fasta("/lustrefs/fadru/projects/asfv-ppi/data/ASFV/Georgia-2007-LR743116_trans-rename.fa")

with open("asfv_genenames.txt", 'w') as f:
    for desc in descs:
        f.write(f"{desc}\n")

