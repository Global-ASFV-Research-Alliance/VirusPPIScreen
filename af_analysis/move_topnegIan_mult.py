import os, glob, pdb
import pandas as pd

neg_PPI_file = "/lustrefs/fadru/projects/asfv-ppi/Output/20240108_Ian_Gold50_paired_unpaired_analysis/Ian_high_score_neg_ctrl.csv"
input_fasta_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Yeast_Gold50_AF-multimer/Gold50-neg_ctrl_short/"
output_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Yeast_Gold50_AF-multimer/high_score_neg_short/"
neg_PPI_df = pd.read_csv(neg_PPI_file)

os.makedirs(output_dir, exist_ok=True)
input_fastas = glob.glob(f"{input_fasta_dir}*.fasta")
for PPI in neg_PPI_df['PPI'].tolist():
    proteinA, proteinB = PPI.split('__')[0], PPI.split('__')[1].split('_')[0]
    for fasta in input_fastas:
        basename = os.path.basename(fasta)
        f_proteinA, f_proteinB = basename.split('__')[0].split('-')[1], basename.split('__')[1].split('.fasta')[0]
        if (proteinA == f_proteinA and proteinB == f_proteinB) or (proteinA == f_proteinB and proteinB == f_proteinA):
            os.system(f"cp {fasta} {output_dir}{basename}") 
