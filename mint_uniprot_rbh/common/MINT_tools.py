import os, sys, pdb, glob
import pandas as pd
import common.REST_api as rest
import common.rbh_tools as rbh 

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

def _read_mint_tsv(mint_db_file, column_names=column_names):
    mint = pd.read_csv(mint_db_file, sep='\t', names=column_names, index_col=None)
    return mint

def filter_mint_sym(mint_db_file, 
            species=['Vaccinia', 'vaccinia', 'Vaccinia virus', 'vaccinia virus'], 
            taxID='10254'):
    """
    This script does the AND search for species OR taxID of the MINT database
    """
    mint = _read_mint_tsv(mint_db_file, column_names=column_names)
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

def filter_mint_asym(mint_db_file, 
            species=['Vaccinia', 'vaccinia', 'Vaccinia virus', 'vaccinia virus'], 
            taxID='10254'):
    """
    This script does the AND search for species OR taxID of the MINT database
    """
    mint = _read_mint_tsv(mint_db_file, column_names=column_names)
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
        if protA or protB:
            filtered_df.loc[index] = row
    return filtered_df

def _crossref_uniprot(mint, 
                      uniprot_proteome_genenames_file,
                      uniprot_fasta_proteome,
                      log_file='',
                      temp='/home/jacob.fenster/tmp/'):
    """
    This cross references UniProtKB accesssion numbers to find common ORF/gene names across species
    It first goes through and sees if the accession numbers and gene names match
    Then it does blastp search of the reference proteome to find top hits from alternative proteomes
    uniprot_proteome_genenames is the previously parsed gene and ORF names of the target proteome file path
    uniprot_fasta_proteome is the fasta formatted UniProt proteome file path
    """
    uniprot_proteome_genenames = pd.read_csv(uniprot_proteome_genenames_file, index_col=0)
    acc_protsA, acc_protsB = [], [] # holds list of accession numbers 
    # get list of protA and protB accession numbers
    for index in mint.index:
        acc_protA, acc_protB = mint.loc[index, 'idA'].split(':')[1], mint.loc[index, 'idB'].split(':')[1]
        acc_protsA.append(acc_protA), acc_protsB.append(acc_protB)
    master_accession = acc_protsA + acc_protsB 
    # search REST api for gene names and ORF names etc. 
    summary = rest.gene_names_from_accession(prot_accessions=master_accession, 
                                             log_file=log_file)
    crossref_mint = mint.copy() 
    for accession in summary.index: 
        try:
            alt_names = summary.loc[accession, 'gene_names'].split(', ')
        except AttributeError:
            alt_names = []
        crossref_list = []
        for WR_accession in uniprot_proteome_genenames.index:
            crossref = 0
            if accession == WR_accession:
                crossref = WR_accession
            if crossref == 0:
                for name in alt_names:
                    try:
                        if name in uniprot_proteome_genenames.loc[WR_accession, 'gene_names']:
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
    # add crossref by blastp
    accessions = []
    for index in crossref_mint.index:
        acc_protA, acc_protB = crossref_mint.loc[index, 'idA'].split(':')[1], crossref_mint.loc[index, 'idB'].split(':')[1]
        if acc_protA == crossref_mint.loc[index, 'crossrefA'] and acc_protB == crossref_mint.loc[index, 'crossrefB']:
            crossref_mint.loc[index, 'genome_match'] = 1
        else: 
            crossref_mint.loc[index, 'genome_match'] = 0
            accessions.append(acc_protA), accessions.append(acc_protB)
    summary = rest.aa_seq_from_accession(prot_accessions=accessions, log_file=log_file)
    os.makedirs(temp, exist_ok=True)
    with open(f'{temp}query.fasta', 'w') as f:
        for accession in summary.index:
            f.write(f">{accession}\n{summary.loc[accession, 'aa_seq']}\n")
    hits = rbh.fbh(f'{temp}query.fasta', uniprot_fasta_proteome, temp,  
                   hit_len_tr=50, qcovs_tr=50, hcovs_tr=50,   
                   outfmt='10 qseqid sseqid evalue pident bitscore qcovs slen qlen length qstart qend sstart send',
                   evalue=0.01)
    # hits is a list of dataframes. Also need to change the vacc proteome file to just have accesssions
    hit_dicts = {}
    for hit in hits:
        if len(hit) > 1:
            # not sure what to do when more than one hit. need to rank them
            breakpoint()
        else:
            h_dict = hit.loc[hit.index[0]].to_dict() # takes first line and converts to dictionary. 
            hit_dicts[h_dict['qseqid']] = h_dict # dictionary with keys of query accession and returns blast dict
    # THis results in a hit dictionary with one hit per accession (key)
    # load back onto crossref_mint
    for index in crossref_mint.index:
        if crossref_mint.loc[index, 'genome_match'] == 0:
            acc_protA, acc_protB = crossref_mint.loc[index, 'idA'].split(':')[1], crossref_mint.loc[index, 'idB'].split(':')[1]
            try:
                crossref_mint.loc[index, 'top_hit_A'] = hit_dicts[acc_protA]['sseqid'].split('|')[1]
                crossref_mint.loc[index, 'pidentA'] = hit_dicts[acc_protA]['pident']
                crossref_mint.loc[index, 'hcovsA'] = hit_dicts[acc_protA]['hcovs']
                crossref_mint.loc[index, 'evalueA'] = hit_dicts[acc_protA]['evalue']
            except KeyError:
                x=1
            try:
                crossref_mint.loc[index, 'top_hit_B'] = hit_dicts[acc_protB]['sseqid'].split('|')[1]
                crossref_mint.loc[index, 'pidentB'] = hit_dicts[acc_protB]['pident']
                crossref_mint.loc[index, 'hcovsB'] = hit_dicts[acc_protB]['hcovs']
                crossref_mint.loc[index, 'evalueB'] = hit_dicts[acc_protB]['evalue']
            except KeyError:
                x=1
    return crossref_mint

def load_mint_to_master(mint, master):
    """
    This takes a mint PPI table (typically from filter and crossref) 
    and converts the PPI to the standard alphabetical format for loading onto a master table
    """
# convert MINT PPIs into standard format 
    ppis = set()
    for index in mint.index:
        if mint.loc[index, 'genome_match'] == 1:
            prot1, prot2 = mint.loc[index, 'crossrefA'], mint.loc[index, 'crossrefB']
        else:
            prot1, prot2 = mint.loc[index, 'top_hit_A'], mint.loc[index, 'top_hit_B']
        if type(prot1) == str and type(prot2) == str:
            protA, protB = min([prot1, prot2]), max([prot1, prot2])
            ppi = f"{protA}__{protB}"
            ppis.add(ppi)
    master['MINT_PPI'] = 0 
    for index in master.index:
        if index in ppis:
            master.loc[index, 'MINT_PPI'] = 1
    return master
