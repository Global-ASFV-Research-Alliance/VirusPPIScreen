# this script is to check to see if all afmult input fastas have been calculated 
import os, sys, glob
from pathlib import Path
import pandas as pd

def _df_not_computed(fasta_in, afmult_out):
    fastas = glob.glob(f"{fasta_in}*.fasta")
    out_folders = glob.glob(f"{afmult_out}*/")
    not_computed = pd.DataFrame(columns=['no_folder_out', 'no_predictions'])
    for fasta in fastas:
        folder_found, predict_found = False, False
        in_ppi = os.path.basename(fasta).split('.')[0]
        for folder in out_folders:
            folder_ppi = Path(folder).name
            if in_ppi == folder_ppi:
                folder_found = True
                if os.path.exists(f"{folder}ranking_debug.json"):
                    predict_found = True
                else:
                    not_computed.loc[in_ppi] = {'no_folder_out': False, 'no_predictions': True}
        if not folder_found: 
            not_computed.loc[in_ppi] = {'no_folder_out': True, 'no_predictions': True}
    return not_computed 


def _delete_empty_folders(not_computed, afmult_out):
    for ppi in not_computed.index:
        os.system(f'rm -r {afmult_out}{ppi}/')
        print(f"deleted folder {afmult_out}{ppi}")

def _copy_not_computed_fastas(not_computed, fasta_in, new_fasta_dir):
    os.makedirs(new_fasta_dir, exist_ok=True)
    for ppi in not_computed.index:
        os.system(f"cp {fasta_in}{ppi}.fasta {new_fasta_dir}")


fasta_in = input(f'Please input the path to the input fasta sequences: ')
afmult_out = input(f'Please input the path to the AFmult output dir: ')
not_computed = _df_not_computed(fasta_in, afmult_out)
print(f"there are {len(not_computed.index)} ppis not computed")
if len(not_computed.index) > 0:
    for ppi in not_computed.index:
        print(f"{ppi} has no predictions and no folder out = {not_computed.loc[ppi, 'no_folder_out']}")
    move_uncomputed = input("Do you want to copy the uncomputed fasta sequences to a new dir? (True/False): ")
    if move_uncomputed:
        new_fasta_dir = input("Please enter the folder path to copy the uncomputed fastas: ")
        _copy_not_computed_fastas(not_computed, fasta_in, new_fasta_dir)
        print(f"copied uncomputed fastas to {new_fasta_dir}")
    delete_empty_folders = input(f"Do you want to delete the folders that have uncomputed models? (True/False): ")
    if delete_empty_folders: 
        _delete_empty_folders(not_computed, afmult_out)
        print(f"deletion of uncomputed model folders complete")
else:
    print('nothing more to do! closing...')
