import os, sys, glob, pdb

# This script is to split the pMSAs into groups so that I can run on multiple nodes
# splitting into 6 tranches

project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/"

tranche1 = f"{project_dir}pmsas/tranche1-pu/"
tranche2 = f"{project_dir}pmsas/tranche2-pu/"
tranche3 = f"{project_dir}pmsas/tranche3-pu/"
tranche4 = f"{project_dir}pmsas/tranche4-pu/"
tranche5 = f"{project_dir}pmsas/tranche5-pu/"
tranche6 = f"{project_dir}pmsas/tranche6-pu/"

os.makedirs(tranche1, exist_ok=True), os.makedirs(tranche2, exist_ok=True), os.makedirs(tranche3, exist_ok=True)
os.makedirs(tranche4, exist_ok=True), os.makedirs(tranche5, exist_ok=True), os.makedirs(tranche6, exist_ok=True)

pmsas = glob.glob(f"{project_dir}pmsas/paired_unpaired/*.a3m")
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

tranche1_pMSAs = glob.glob(f"{tranche1}*.a3m")
tranche2_pMSAs = glob.glob(f"{tranche2}*.a3m")
tranche3_pMSAs = glob.glob(f"{tranche3}*.a3m")
tranche4_pMSAs = glob.glob(f"{tranche4}*.a3m")
tranche5_pMSAs = glob.glob(f"{tranche5}*.a3m")
tranche6_pMSAs = glob.glob(f"{tranche6}*.a3m")
a = set()
for msa in tranche1_pMSAs:
    a.add(msa)
for msa in tranche2_pMSAs:
    a.add(msa)
for msa in tranche3_pMSAs:
    a.add(msa)

for msa in tranche4_pMSAs:
    a.add(msa)

for msa in tranche5_pMSAs:
    a.add(msa)
for msa in tranche6_pMSAs:
    a.add(msa)
breakpoint()
