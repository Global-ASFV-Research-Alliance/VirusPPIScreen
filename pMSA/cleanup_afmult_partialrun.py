# this script is to cleanup folders from an aborted afmult run 
# that uses the AF_multimer_msa.py script to move msas from run
import os, glob, pdb

def cleanup_afmult_output(output_dir):
    runs = glob.glob(f"{output_dir}*/")
    for run in runs:
        if os.path.exists(f"{run}ranking_debug.json"):
            continue
        else:
            print(f"rm -r {run}")
            os.system(f"rm -r {run}")

def cleanup_input_fasta(input_dir, move_input_dir, output_dir):
    os.makedirs(move_input_dir, exist_ok=True)
    input_fastas = glob.glob(f"{input_dir}*.fasta")
    output_models = glob.glob(f"{output_dir}*/")
    output_model_names = []
    for model in output_models:
        output_model_names.append(os.path.basename(model.rstrip(os.sep)))
    for fasta in input_fastas:
        name = os.path.basename(fasta).split('.')[0]
        if name in output_model_names:
            print(f"mv {fasta} {move_input_dir}{name}.fasta")
            os.system(f"mv {fasta} {move_input_dir}{name}.fasta")

input_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/alphafold_multimer/high_score_heterodimer/"
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/ASFV/alphafold-multimer/high_score_heterodimer/"
move_input_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/alphafold_multimer/high_score_heterodimer/Already_Computed/"

cleanup_afmult_output(output_dir)
#cleanup_input_fasta(input_dir, move_input_dir, output_dir)
