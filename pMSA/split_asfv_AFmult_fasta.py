import os, sys, glob, pdb

# This script is to split the pMSAs into groups so that I can run on multiple nodes
# splitting into 6 tranches

project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/alphafold_multimer/mid_high_score_heterodimer/"

tranche1 = f"{project_dir}tranche1/"
tranche2 = f"{project_dir}tranche2/"
tranche3 = f"{project_dir}tranche3/"
tranche4 = f"{project_dir}tranche4/"
tranche5 = f"{project_dir}tranche5/"
tranche6 = f"{project_dir}tranche6/"

os.makedirs(tranche1, exist_ok=True), os.makedirs(tranche2, exist_ok=True), os.makedirs(tranche3, exist_ok=True)
os.makedirs(tranche4, exist_ok=True), os.makedirs(tranche5, exist_ok=True), os.makedirs(tranche6, exist_ok=True)

pmsas = glob.glob(f"{project_dir}*.fasta")
num_pmsas = len(pmsas)
print(f"there are {num_pmsas} pMSAs. Splitting into six tranches")
breakpoint()

if True:
    for i in range(0, int( num_pmsas/6)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche1}{basename}")
    for i in range(int(num_pmsas/6), int(num_pmsas/6 * 2)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche2}{basename}")
    for i in range(int(num_pmsas/6 * 2), int(num_pmsas/6 * 3)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche3}{basename}")
    for i in range(int(num_pmsas/6 * 3), int(num_pmsas/6 * 4)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche4}{basename}")
    for i in range(int(num_pmsas/6 * 4), int(num_pmsas/6 * 5)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche5}{basename}")
    for i in range(int(num_pmsas/6 * 5),int( num_pmsas)):
        basename = os.path.basename(pmsas[i])
        os.system(f"cp {pmsas[i]} {tranche6}{basename}")

