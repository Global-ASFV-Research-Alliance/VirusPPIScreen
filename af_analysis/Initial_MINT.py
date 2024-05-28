# this was moved to mint uniprot rbh
# This script is to calculate precision recall curves for the vaccinia dataset 
import os, sys, pdb, glob
import analysis.AFppi as ppi
import pandas as pd
import numpy as np
sys.path.append("/home/jacob.fenster/scripts/api/common/")
import REST_api as rest
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar

def _filter_sym(mint, 
            species=['Vaccinia', 'vaccinia', 'Vaccinia virus', 'vaccinia virus'], 
            taxID='10254'):
    """
    This script does the AND search for species OR taxID
    """
    filtered_df = pd.DataFrame(columns=mint.columns)
    for index in mint.index:
        row = mint.iloc[index].to_dict()
        protA, protB = False, False # consensus switch to add to filtered df
        for name in species:
            if name in row['taxidA']:
                protA = True
            if name in row['taxidB']:
                protB = True
        if taxID in row['taxidA']:
            protA = True
        if taxID in row['taxidB']:
            protB = True
        if protA and protB:
            filtered_df.loc[index] = row
    return filtered_df

def _crossref_uniprot(mint, uniprot_proteome, api=True, summary=None):
    """
    This cross references UniProtKB accesssion numbers to find common ORF/gene names across species
    """
    log_file=f"{vaccinia_out}MINT_uniprot_crossref.log"
    acc_protsA, acc_protsB = [], [] # holds list of accession numbers 
    # get list of protA and protB accession numbers
    for index in mint.index:
        acc_protA, acc_protB = mint.loc[index, 'idA'].split(':')[1], mint.loc[index, 'idB'].split(':')[1]
        acc_protsA.append(acc_protA), acc_protsB.append(acc_protB)
    master_accession = acc_protsA + acc_protsB 
    if api:
        summary = rest.gene_names_from_accession(prot_accessions=master_accession, 
                                                 log_file=log_file)
    else:
        summary = pd.read_csv(summary, index_col=0)
    crossref_mint = mint.copy() 
    for accession in summary.index: 
        try:
            alt_names = summary.loc[accession, 'gene_names'].split(', ')
        except AttributeError:
            alt_names = []
        crossref_list = []
        for WR_accession in uniprot_proteome.index:
            crossref = 0
            if accession == WR_accession:
                crossref = WR_accession
            if crossref == 0:
                for name in alt_names:
                    try:
                        if name in uniprot_proteome.loc[WR_accession, 'gene_names']:
                            crossref = WR_accession
                    except TypeError:
                        os.system(f'echo "{WR_accession} does not have crossref" >> {log_file}')
            if crossref != 0:
                crossref_list.append(crossref)
        if len(crossref_list) == 0:
            os.system(f'echo "{accession} has no crossreference" >> {log_file}')
        elif len(crossref_list) > 1:
            os.system(f'echo "{accession} has multiple crossreferences: {crossref_list}" >> {log_file}')
        else: 
            summary.loc[accession, 'crossref'] = crossref_list[0]
    # now input crossref in MINT df
    for index in crossref_mint.index:
        acc_protA, acc_protB = mint.loc[index, 'idA'].split(':')[1], mint.loc[index, 'idB'].split(':')[1]
        crossref_A, crossref_B = summary.loc[acc_protA, 'crossref'], summary.loc[acc_protB, 'crossref']
        crossref_mint.loc[index, 'crossrefA'], crossref_mint.loc[index, 'crossrefB'] = crossref_A, crossref_B
    return crossref_mint

# column names for MINT database
column_names = [
        'idA', # identifier A 
        'idB', # identifier B
        'altIDsA', # Alternative Identifiers A
        'altIDsB', # Alterntaive identifiers B
        'aliasA', # alias A
        'aliasB', # alias B
        'detmethod', # iteraction detection method
        'pubauth', # publication first author
        'pubid', # publication identifier
        'taxidA', # Tax ID interactor A: Species or tax ID
        'taxidB', # Tax ID interactor B: Species or tax ID
        'type', # interaction type(s)
        'source', # Interaction source(s)
        'interaction_id', # interaction identifier
        'confidence'
        ]
vaccinia_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
descs, seqs = bpar.read_fasta("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta")
vacc_proteome = {}
for desc, seq in zip(descs, seqs):
    accession = desc.split('|')[1]
    vacc_proteome[accession] = seq

if False:
    mint = pd.read_csv("/lustrefs/fadru/projects/asfv-ppi/data/MINT/MINT_db.tsv", sep='\t', names=column_names, index_col=None)
# filter for vaccinia-vaccinia virus PPI
    if False:
        filtered_mint = _filter_sym(mint)
        filtered_mint.to_csv(f"{vaccinia_out}MINT_vaccinia-vaccinia_ppi.csv")
    filtered_mint = pd.read_csv(f"{vaccinia_out}MINT_vaccinia-vaccinia_ppi.csv", index_col=0)
# cross reference with UniProtKB to get Vaccinia-WR specific PPI
    uniprot_vaccWR = pd.read_csv(f"/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_virus-WR_UniProt_genenames.csv", index_col=0)
    crossref_mint = _crossref_uniprot(filtered_mint, uniprot_vaccWR,
                                      api=False, summary=f"{vaccinia_out}MINT_vaccinia_genenames.csv")      
    crossref_mint.to_csv(f"{vaccinia_out}MINT_vacciniaWR-vacciniaWR_crossref.csv")
# convert MINT PPIs into standard format 
    ppis = set()
    for index in crossref_mint.index:
        prot1, prot2 = crossref_mint.loc[index, 'crossrefA'], crossref_mint.loc[index, 'crossrefB']
        if type(prot1) == str and type(prot2) == str:
            protA, protB = min([prot1, prot2]), max([prot1, prot2])
            ppi = f"{protA}__{protB}"
            ppis.add(ppi)
    vacc_master = pd.read_csv(f"{vaccinia_out}vacc_master_partial8.csv", index_col=0)
    vacc_master['MINT_PPI'] = 0 
    for index in vacc_master.index:
        if index in ppis:
            vacc_master.loc[index, 'MINT_PPI'] = 1
    vacc_master.to_csv(f"{vaccinia_out}vacc_master_partial9.csv")
# make AFmult fastas for the new PPI gold set 

# filter for in MINT gold std and NOT AFmult calc
    afmult_fasta_out = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/MINT_positives/"
    vacc_afmult_run = vacc_master[(vacc_master['MINT_PPI'] == 1) & (vacc_master['iptm+ptm'].isna())]
    for ppi in vacc_afmult_run.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1] 
        with open(f"{afmult_fasta_out}{ppi}.fasta", 'w') as f:
            f.write(f">{protA}\n{vacc_proteome[protA]}\n>{protB}\n{vacc_proteome[protB]}")
# Now create preformance with genomewide information
vacc_master = pd.read_csv(f"{vaccinia_out}vacc_master_partial9.csv", index_col=0)
# create positive control set
vacc_positive = vacc_master[vacc_master['MINT_PPI'] == 1]
# create negative control set 
vacc_neg_ppi = vacc_master[vacc_master['MINT_PPI'] == 0]

threshold_arr = np.arange(0, 1, .01)
column = 'Humphreys_maxcontact'
precision_list, recall_list, df = ppi.precision_recall_curve(vacc_positive, vacc_neg_ppi, threshold_arr, column)
breakpoint()
# now calculate HAF + AFmult preformance. Assumes that all uncalculated AFmult are 0
vacc_above_tr = vacc_master[vacc_master[column] >= 0.8]
# create fasta files for afmult of nan values
if False:
    vacc_above_tr_nan = vacc_above_tr[vacc_above_tr['iptm+ptm'].isna()]
    afmult_fasta_out = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/afmult_ian-tophits/leftovers/"
    for ppi in vacc_above_tr_nan.index:
        protA, protB = ppi.split('__')[0], ppi.split('__')[1]
        with open(f"{afmult_fasta_out}{ppi}.fasta", 'w') as f:
            f.write(f">{protA}\n{vacc_proteome[protA]}\n>{protB}\n{vacc_proteome[protB]}")
vacc_afmult_curve = vacc_master.copy()
vacc_afmult_curve = vacc_afmult_curve.where(
    # https://note.nkmk.me/en/python-pandas-where-mask/
