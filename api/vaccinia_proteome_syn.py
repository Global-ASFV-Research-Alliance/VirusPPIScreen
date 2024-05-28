import os, pdb
import pandas as pd
import common.REST_api as rest
import common.bio_parsers as bpar 

output_dir = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/"
fasta_output_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/AFmult/"
vaccinia_uniprot_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
mccraith2000_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000_vaccinia_y2h.csv"

proteome_df = bpar.parse_uniprot_proteome(vaccinia_uniprot_proteome)
if False: # pull from UniProt REST
    prot_accessions = proteome_df.index.tolist()
    vacc_names = rest.gene_names_from_accession(prot_accessions)
    vacc_names.to_csv(f"{output_dir}vaccinia_virus-WR_UniProt_genenames.csv")
mccraith_df = pd.read_csv(mccraith2000_file)
vacc_names = pd.read_csv(f"{output_dir}vaccinia_virus-WR_UniProt_genenames.csv", index_col=0)

# generate AFmult PPIs from the McCraith 2000 ppi data
os.system(f"> {fasta_output_dir}failed_fastas.log")
for index in mccraith_df.index:
    protA, protB = mccraith_df.loc[index, 'protA'], mccraith_df.loc[index, 'protB']
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
        with open(f"{fasta_output_dir}{protA}__{protB}.fasta", 'w') as f:
            if protA == protB:
                f.write(f">{protA}__A\n{proteome_df.loc[protA_accession[0], 'sequence']}\n")
                f.write(f">{protB}__B\n{proteome_df.loc[protB_accession[0], 'sequence']}")
            else:
                f.write(f">{protA}\n{proteome_df.loc[protA_accession[0], 'sequence']}\n")
                f.write(f">{protB}\n{proteome_df.loc[protB_accession[0], 'sequence']}")
    else:
        os.system(f'echo "PPI {protA}__{protB} failed. {protA} has {len(protA_accession)} accesssion(s). {protB} has {len(protB_accession)} accessions(s)." >> {fasta_output_dir}failed_fastas.log')
os.system(f'echo "fasta generation complete" >> {fasta_output_dir}failed_fastas.log')
