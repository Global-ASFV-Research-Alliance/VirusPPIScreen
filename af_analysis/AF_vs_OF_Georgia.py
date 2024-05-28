import os, glob, pdb
import pandas as pd
import analysis.AFppi as ppi

# this script calculates the average pLDDT of the AF and OF calculated models of corresponding proteins 
af_asfv_georiga_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_monomer_georgia/"
of_asfv_georiga_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/OF_ASFV_monomers/"
af_paths = glob.glob(f"{af_asfv_georiga_dir}*.pdb")
of_paths = glob.glob(f"{of_asfv_georiga_dir}*.pdb")

output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/20240117_AF_vs_OF_ASFV_monomer/"

summary = pd.DataFrame(columns=['af-avg_plddt', 'of-avg_plddt'])
# parse OF ASFV georiga monomers
for file in of_paths:
    basename = os.path.basename(file).split('.pdb')[0]
    protein = basename.split('Georgia-2007_')[1]
    summary.loc[protein, 'of-avg_plddt'] = ppi.avg_plddt_from_pdb(file)
for file in af_paths:
    basename = os.path.basename(file).split('.pdb')[0]
    protein = basename.replace(' ', '_')
    summary.loc[protein, 'af-avg_plddt'] = ppi.avg_plddt_from_pdb(file)
breakpoint()
summary['of_minus_af'] = summary['of-avg_plddt'] - summary['af-avg_plddt']

summary.to_csv(f"{output_dir}AF_vs_OF_average_plddt_data.csv")

