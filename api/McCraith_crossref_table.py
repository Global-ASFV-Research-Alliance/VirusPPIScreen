# This script is to generate a table of McCraith 2000 gold standards for publication
import os, pdb
import pandas as pd
import numpy as np
import common.bio_parsers as bpar

output_dir = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/"
vacc_names = pd.read_csv(f"{output_dir}vaccinia_virus-WR_UniProt_genenames.csv", index_col=0)
fasta_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/McCraith2000/"
vaccinia_uniprot_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
mccraith2000_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000_vaccinia_y2h.csv"

# make proteome dict
vaccinia_wr_dict = {}
descs, seqs = bpar.read_fasta(vaccinia_uniprot_proteome)
for desc, seq in zip(descs, seqs):
    accession = desc.split('|')[1]
    vaccinia_wr_dict[accession] = seq

table = pd.DataFrame(columns=['protA_desc', 'protB_desc', 'protA_seq', 'protB_seq'])
mccraith_df = pd.read_csv(mccraith2000_file)
for index in mccraith_df.index:
    protA, protB = mccraith_df.loc[index, 'protA'], mccraith_df.loc[index, 'protB']
    ppi = f"{protA}__{protB}"
    protA_accession, protB_accession = [], []
    for accession in vacc_names.index:
        try:
            for name in vacc_names.loc[accession, 'gene_names'].split(', '):
                if protA == name:
                    protA_accession.append(accession)
                if protB == name: 
                    protB_accession.append(accession)
        except TypeError:
            continue
        except AttributeError:
            continue
    if len(protA_accession) == 1 and len(protB_accession) == 1:
        protA_seq, protB_seq = vaccinia_wr_dict[protA_accession[0]], vaccinia_wr_dict[protB_accession[0]] 
        table.loc[ppi] = {'protA_desc': mccraith_df.loc[index, 'protA_desc'],
                          'protB_desc': mccraith_df.loc[index, 'protB_desc'],
                          'protA_seq': protA_seq, 
                          'protB_seq': protB_seq}
    elif len(protA_accession) == 1:
        protA_seq, protB_seq = vaccinia_wr_dict[protA_accession[0]], np.nan 
        table.loc[ppi] = {'protA_desc': mccraith_df.loc[index, 'protA_desc'],
                          'protB_desc': mccraith_df.loc[index, 'protB_desc'],
                          'protA_seq': protA_seq, 
                          'protB_seq': protB_seq}
    elif len(protB_accession) == 1:
        protA_seq, protB_seq = np.nan, vaccinia_wr_dict[protB_accession[0]]
        table.loc[ppi] = {'protA_desc': mccraith_df.loc[index, 'protA_desc'],
                          'protB_desc': mccraith_df.loc[index, 'protB_desc'],
                          'protA_seq': protA_seq, 
                          'protB_seq': protB_seq}
    else: 
        table.loc[ppi] = {'protA_desc': mccraith_df.loc[index, 'protA_desc'],
                          'protB_desc': mccraith_df.loc[index, 'protB_desc']}
        
breakpoint() 


