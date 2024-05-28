# This script is to cross reference each viral genome in the NCBI database with its taxID and other metadata
import os, sys, glob, pdb, time
import pandas as pd
import numpy as np
import common.Entrez_NCBI as entrez

# data input
api_key = "46f09d13b8f485ab42dabf140b3a09dd9c09"
api_key_usda = '332aca54a1bc5882797d51672572237b6c08'
query_all_NCBI = 'Viruses[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[prop]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[Filter])'
query_all_virus = 'Viruses[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[PROP]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[FILTER]+OR+complete+genome[TITL])'
query_ASFV = 'Viruses[Organism]+AND+African+swine+fever+virus[Organism]+NOT+cellular+organisms[ORGN]+NOT+wgs[PROP]+NOT+gbdiv+syn[PROP]+AND+(srcdb_refseq[PROP]+OR+nuccore+genome+samespecies[FILTER]+OR+complete+genome[TITL])'
db = 'nuccore'
exp_tag = '20240411_all_virus'
output_dir = '/home/jacob.fenster/scripts/api/data/'

# esearch nuccore and return/save list of accession numbers
if False:
    search_Ids = entrez.esearch_id_list(query_all_virus, db, api_key, exp_tag, output_dir, retmax=10000)
# load prevsiouly queried search_Ids (accession numbers)
search_Ids_df = pd.read_csv('/home/jacob.fenster/scripts/api/data/20240411_all_virus_searchIds.csv')
search_Ids = search_Ids_df['acession'].tolist()
# pull all taxonomy lineage data from list of viral genomes
if False:
    tax_master_list = []
    tax_master_list = entrez.elink_nuccore_taxonomy(search_Ids, tax_master_list, api_key, exp_tag, output_dir, begin=0, retmax=200)
    df_columns = ['ParentTaxId', 'superkingdom', 'clade', 'kingdom', 'phylum', 'subphylum', 'class', 'order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'no rank', 'entry_rank']
    tax_df = pd.DataFrame(columns=df_columns)
    tax_df.index.name = 'TaxId'
    tax_df = entrez.efetch_taxonomy(tax_master_list, tax_df, api_key, output_dir, exp_tag, begin=0, retmax=200)
    tax_df.to_csv(f"{output_dir}/{exp_tag}_taxonomy.csv")
# tax_df has each virus taxon with lineage information. Previously computed, code above
tax_df = pd.read_csv(f"/home/jacob.fenster/scripts/api/data/20231106_all_virus_taxonomy.csv", index_col=0)
# try with subset
search_Ids = search_Ids[:200]
breakpoint()
timei = time.time()
entrez.elink_single_nuccore_taxonomy(search_Ids, api_key, exp_tag, output_dir, begin=0, end=len(search_Ids), savepoint=5, log=f'{output_dir}{exp_tag}_elink.log')
timef = time.time()
print(f'time elapsed for {len(search_Ids)} is {timef-timei:.1f}')
breakpoint()




# pull all metadata from nuccore databse of each viral genome
if False:
    nuccore_meta_df = pd.DataFrame()
    nuccore_meta_df.index.name = 'acession'
    nuccore_meta_df = efetch_metatadata_nuccore(nuccore_meta_df, search_Ids, api_key, output_dir, exp_tag, begin=0, retmax=200)
    nuccore_meta_df.to_csv(f"{output_dir}/{exp_tag}_nuccore_meta.csv")
    pdb.set_trace()
