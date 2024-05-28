import os, pdb
import pandas as pd
import analysis.plot as jplot

output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20240108_Ian_Gold50_paired_unpaired_analysis/"
hineg_violin = pd.read_csv(f"{output_dir}hiscore_neg_violin_data.csv", index_col=0)
breakpoint()
jplot.stripplot(hineg_violin, x='experiment', y='score', hue=None, output_dir=output_dir, tag='high_score_neg_stripplot')
