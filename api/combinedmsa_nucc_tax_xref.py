# This script is to cross reference the viral genomes in the ASFV/vaccinia combined MSAs to the tax db 
#
import os, sys, glob, pdb, time
import pandas as pd
import numpy as np
import common.Entrez_NCBI as entrez
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar

# data input
api_key = "46f09d13b8f485ab42dabf140b3a09dd9c09"
api_key_usda = '332aca54a1bc5882797d51672572237b6c08'
db = 'nuccore'
output_dir = '/home/jacob.fenster/scripts/api/data/msa_only/'

### Switches
elink_msas = False

### extract all db accessions from the combined MSAs and save to csv
if elink_msas:
# ASFV
    asfv_combined_msas = glob.glob("/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/combined_msas/*.a3m")

    asfv_accessions = set()
    for msa in asfv_combined_msas:
        descs, seqs = bpar.read_fasta(msa)
        for desc in descs:
            if desc.startswith('lcl|'):
                accession = desc.split('lcl|')[1].split('_prot')[0]
                asfv_accessions.add(accession)
    asfv_accessions = list(asfv_accessions)
    with open(f'{output_dir}asfv_accession_list.csv', 'w') as f:
        f.write('accession,\n')
        for asfv_accession in asfv_accessions:
            f.write(f'{asfv_accession},\n')
# Vacc
    vacc_combined_msas = glob.glob("/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFIan2021/vaccinia-vaccinia-mega/combined_msas/*.a3m")
    vacc_accessions = set()
    for msa in vacc_combined_msas:
        descs, seqs = bpar.read_fasta(msa)
        for desc in descs:
            if desc.startswith('lcl|'):
                accession = desc.split('lcl|')[1].split('_prot')[0]
                vacc_accessions.add(accession)
    vacc_accessions = list(vacc_accessions)
    with open(f'{output_dir}vacc_accession_list.csv', 'w') as f:
        f.write('accession,\n')
        for vacc_accession in vacc_accessions:
            f.write(f'{vacc_accession},\n')
# ASFV xref
    exp_tag = 'asfv_combined_msas'
    begin, end = 0, len(asfv_accessions)
    asfv_xref = entrez.elink_single_nuccore_taxonomy(asfv_accessions, api_key, exp_tag, output_dir, 
                                                     begin=begin, end=end, savepoint=6000, log=f'{output_dir}{exp_tag}_elink.log')
# Vacc xref
vacc_accessions = pd.read_csv(f'{output_dir}vacc_accession_list.csv', index_col=0)
exp_tag = 'vacc_combined_msas'
begin, end = 0, len(vacc_accessions)
vacc_xref = entrez.elink_single_nuccore_taxonomy(vacc_accessions, api_key, exp_tag, output_dir, 
                                                 begin=begin, end=end, savepoint=6000, log=f'{output_dir}{exp_tag}_elink.log')



## extract errors 

def _extract_errors(log_file):
    # incomplete 
    with open(log_file, 'r') as f:
        lines = f.readlines()
        breakpoint()
        accession = ''
        for i in range(0, len(lines)):
            if lines[i].startswith('ERROR multiple'):
                breakpoint()
                xml = ''
                accession = lines[i].split('xrefs: ')[1].split('.')[0]
                j = i + 1 
                while not lines[j].strip().startswith('<') or not lines.strip().startswith(''):
                    xml += lines[j]
                    j += 1
                xref = root.find('./Link')


asfv_log = '/home/jacob.fenster/scripts/api/data/msa_only/asfv_combined_msas_elink.log'
vacc_log = '/home/jacob.fenster/scripts/api/data/msa_only/vacc_combined_msas_elink.log'
