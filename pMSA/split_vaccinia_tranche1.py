import os, sys, glob, pdb

project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFIan2021/vaccinia-vaccinia-mega/"

original = f"{project_dir}pmsas/tranche1-pu/"
output_og = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/Ian2021/homo_heterodimer/tranche1/"

tranche1 = f"{project_dir}pmsas/tranche1_1-pu/"
tranche2 = f"{project_dir}pmsas/tranche1_2-pu/"
tranche3 = f"{project_dir}pmsas/tranche1_3-pu/"
tranche4 = f"{project_dir}pmsas/tranche1_4-pu/"
tranche5 = f"{project_dir}pmsas/tranche1_5-pu/"
tranche6 = f"{project_dir}pmsas/tranche1_6-pu/"

os.makedirs(tranche1, exist_ok=True), os.makedirs(tranche2, exist_ok=True), os.makedirs(tranche3, exist_ok=True)
os.makedirs(tranche4, exist_ok=True), os.makedirs(tranche5, exist_ok=True), os.makedirs(tranche6, exist_ok=True)


all_pmsas = glob.glob(f"{original}*.a3m")
output_files = glob.glob(f"{output_og}*.npz")
pmsas = []
breakpoint()
for pmsa in all_pmsas:
    basename = os.path.basename(pmsa).split('.')[0] 
    name = f"{basename.split('-AA')[0]}__{basename.split('__')[1].split('-AA')[0]}"
    output_name = f"{output_og}{name}_info_model_3_ptm.npz"
    if output_name not in output_files:
        pmsas.append(pmsa)

num_pmsas = len(pmsas)
breakpoint()

 
if True:
    for i in range(0, int(num_pmsas/6)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche1}{basename}")
    for i in range(int(num_pmsas/6), int(num_pmsas/6 * 2)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche2}{basename}")
    for i in range(int(num_pmsas/6 * 2), int(num_pmsas/6 * 3)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche3}{basename}")
    for i in range(int(num_pmsas/6 * 3), int(num_pmsas/6 * 4)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche4}{basename}")
    for i in range(int(num_pmsas/6 * 4), int(num_pmsas/6 * 5)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche5}{basename}")
    for i in range(int(num_pmsas/6 * 5),int( num_pmsas)):
        basename = os.path.basename(pmsas[i])
        chain1, chain2  = int(basename.split('__')[0].split('-AA')[1]), int(basename.split('__')[1].split('-AA')[1].split('-')[0])
        if chain1 > 1200 or chain2 > 1200:
            os.system(f"cp {pmsas[i]} {large_peptides}{basename}")
        else:
            os.system(f"cp {pmsas[i]} {tranche6}{basename}")
