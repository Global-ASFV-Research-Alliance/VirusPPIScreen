# this script is to parse the file names of the ASFV-ASFV all vs all hetero-homodimer
# set for running on Ian 2021 AF script
import sys, os

file_path = sys.argv[1]
basename = os.path.basename(file_path).split('.')[0]

gene1 = basename.split('-AA')[0] 
gene2 = basename.split('__')[1].split('-AA')[0]
chain1 = basename.split('__')[0].split('-AA')[1]
chain2 = basename.split('__')[1].split('-AA')[1].split('-')[0]

print(f"gene1={gene1}")
print(f"gene2={gene2}")
print(f"chain1={chain1}")
print(f"chain2={chain2}")
