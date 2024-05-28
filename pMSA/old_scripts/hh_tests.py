import os, pdb
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa

file = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/Ian2021/tmp/aln/A22R-unaligned.fasta"
descs, seqs = bpar.read_fasta(file)
for desc, seq in zip(descs, seqs):
    if '.' in seq:
        breakpoint()
    if '-' in seq:
        breakpoint()
    if 'X' in seq:
        breakpoint()
    if 'J' in seq:
        breakpoint()
    if 'B' in seq:
        breakpoint()
    if 'Z' in seq:
        breakpoint()
