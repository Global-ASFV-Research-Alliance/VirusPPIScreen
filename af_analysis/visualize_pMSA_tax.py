import os, sys, glob
import pandas as pd
import numpy as np
import analysis.plot as plot
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar

def _xref_tax_tree(accessions, xref, tax_tree):
    # accessions xref and tax_tree are both df
    xref_full = pd.DataFrame(columns=xref.columns.tolist()+tax_tree.columns.tolist())
    xref_full = xref_full.astype({'taxID': 'Int64', 'ParentTaxId': 'Int64'})
    for accession in accessions.index:
        row = {}
        # add taxID xref
        try:
            taxid = int(xref.loc[accession, 'taxID'])
            row['taxID'] = taxid
        except KeyError:
            # set to nan if no xref and continue 
            xref_full.loc[accession] = np.nan
            continue
        except ValueError:
            if ';' in xref.loc[accession, 'taxID']:
                taxid = int(xref.loc[accession, 'taxID'].split(';')[0])
                row['taxID'] = taxid
            else:
                breakpoint()
        # add tax tree xref
        try:
            tree_row = tax_tree.loc[taxid].to_dict()
            row = row | tree_row
            xref_full.loc[accession] = row
        except KeyError:
            xref_full.loc[accession] = row
    xref_full = xref_full.astype({'taxID': 'Int64', 'ParentTaxId': 'Int64'})
    return xref_full

def _process_column_piechart(df, column):
    # fill np.nan values
    df = df.fillna(value={column: 'other'})
    # get value counts 
    counts = df[column].value_counts()
    labels, data = counts.index.tolist(), counts.tolist()
    return labels, data 


def _extract_accession_pmsa(pmsa_dir, output_csv):
    pmsas = glob.glob(f'{pmsa_dir}*.a3m')
    accessions = set()
    for pmsa in pmsas:
        descs, seqs = bpar.read_fasta(pmsa)
        for desc in descs:
            try: 
                if len(desc.split('.')) == 3:
                    accession = desc.split('.')[0]
                    accessions.add(accession)
            except:
                print(f'ERROR: {desc}')
                continue 
    with open(output_csv, 'w') as f:
        f.write('accession,\n')
        for accession in accessions:
            f.write(f'{accession},\n')
    accession_df = pd.read_csv(output_csv, index_col=0)
    return accession_df

def _remove_accession_version(df):
    for index in df.index:
        new_index = index.split('.')[0]
        df = df.rename(index={index: new_index})
    return df
    
asfv_out = '/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/'
vacc_out = '/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/'
asfv_xref = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/asfv_combined_msas_taxId_xref_final.csv', index_col=0)
vacc_xref = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/vacc_combined_msas_taxId_xref_final.csv', index_col=0)
viral_tax_tree = pd.read_csv('/lustrefs/fadru/projects/asfv-ppi/data/NCBI/20231106_all_virus_taxonomy.csv', index_col=0)

if False:
    viral_tax_tree = pd.read_csv('/lustrefs/fadru/projects/asfv-ppi/data/NCBI/20231106_all_virus_taxonomy.csv', index_col=0)

    asfv_accessions = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/asfv_accession_list.csv', index_col=0)
    asfv_xref = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/asfv_combined_msas_taxId_xref_final.csv', index_col=0)

    vacc_xref = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/vacc_combined_msas_taxId_xref_final.csv', index_col=0)
    vacc_accessions = pd.read_csv('/home/jacob.fenster/scripts/api/data/msa_only/vacc_accession_list.csv', index_col=0)

    asfv_xref = _xref_tax_tree(asfv_accessions, asfv_xref, viral_tax_tree)
    vacc_xref = _xref_tax_tree(vacc_accessions, vacc_xref, viral_tax_tree)
    asfv_xref.to_csv(f'{asfv_out}asfv_tax_xref_combined_msas.csv')
    vacc_xref.to_csv(f'{vacc_out}vacc_tax_xref_combined_msas.csv')

    asfv_xref = pd.read_csv(f'{asfv_out}asfv_tax_xref_combined_msas.csv', index_col=0)
    vacc_xref = pd.read_csv(f'{vacc_out}vacc_tax_xref_combined_msas.csv', index_col=0)
### process for pie charts
    labels, data = _process_column_piechart(asfv_xref, 'kingdom')
    kingdom_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    kingdom_df.to_csv(f'{asfv_out}asfv_kingdom_msa_tax_data.csv')
    labels, data = _process_column_piechart(asfv_xref, 'phylum')
    phylum_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    phylum_df.to_csv(f'{asfv_out}asfv_phylum_msa_tax_data.csv')
    labels, data = _process_column_piechart(asfv_xref, 'class')
    class_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    class_df.to_csv(f'{asfv_out}asfv_class_msa_tax_data.csv')
    labels, data = _process_column_piechart(asfv_xref, 'species')
    species_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    species_df.to_csv(f'{asfv_out}asfv_species_msa_tax_data.csv')


    labels, data = _process_column_piechart(vacc_xref, 'kingdom')
    kingdom_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    kingdom_df.to_csv(f'{vacc_out}vacc_kingdom_msa_tax_data.csv')
    labels, data = _process_column_piechart(vacc_xref, 'phylum')
    phylum_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    phylum_df.to_csv(f'{vacc_out}vacc_phylum_msa_tax_data.csv')
    labels, data = _process_column_piechart(vacc_xref, 'class')
    class_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    class_df.to_csv(f'{vacc_out}vacc_class_msa_tax_data.csv')
    labels, data = _process_column_piechart(vacc_xref, 'species')
    species_df = pd.DataFrame(data={'labels':labels, 'counts': data})
    species_df.to_csv(f'{vacc_out}vacc_species_msa_tax_data.csv')

### extract paired tax info
# ASFV
asfv_pmsas = '/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/pmsas/paired/'
asfv_paired_accessions = _extract_accession_pmsa(asfv_pmsas, f'{asfv_out}asfv_paired_accession_list.csv')

vacc_pmsas = '/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFIan2021/vaccinia-vaccinia-mega/pmsas/paired/'
vacc_paired_accessions =  _extract_accession_pmsa(vacc_pmsas, f'{vacc_out}vacc_paired_accession_list.csv')
asfv_xref = _remove_accession_version(asfv_xref)
vacc_xref = _remove_accession_version(vacc_xref)
breakpoint()

asfv_paired_xref = _xref_tax_tree(asfv_paired_accessions, asfv_xref, viral_tax_tree)
vacc_paired_xref = _xref_tax_tree(vacc_paired_accessions, vacc_xref, viral_tax_tree)

asfv_paired_xref.to_csv(f'{asfv_out}asfv_paired_tax_xref_combined_msas.csv')
vacc_paired_xref.to_csv(f'{vacc_out}vacc_paired_tax_xref_combined_msas.csv')

labels, data = _process_column_piechart(asfv_paired_xref, 'kingdom')
kingdom_df = pd.DataFrame(data={'labels':labels, 'counts': data})
kingdom_df.to_csv(f'{asfv_out}asfv_paired_kingdom_msa_tax_data.csv')
labels, data = _process_column_piechart(asfv_paired_xref, 'phylum')
phylum_df = pd.DataFrame(data={'labels':labels, 'counts': data})
phylum_df.to_csv(f'{asfv_out}asfv_paired_phylum_msa_tax_data.csv')
labels, data = _process_column_piechart(asfv_paired_xref, 'class')
class_df = pd.DataFrame(data={'labels':labels, 'counts': data})
class_df.to_csv(f'{asfv_out}asfv_paired_class_msa_tax_data.csv')
labels, data = _process_column_piechart(asfv_paired_xref, 'species')
species_df = pd.DataFrame(data={'labels':labels, 'counts': data})
species_df.to_csv(f'{asfv_out}asfv_paired_species_msa_tax_data.csv')


labels, data = _process_column_piechart(vacc_paired_xref, 'kingdom')
kingdom_df = pd.DataFrame(data={'labels':labels, 'counts': data})
kingdom_df.to_csv(f'{vacc_out}vacc_paired_kingdom_msa_tax_data.csv')
labels, data = _process_column_piechart(vacc_paired_xref, 'phylum')
phylum_df = pd.DataFrame(data={'labels':labels, 'counts': data})
phylum_df.to_csv(f'{vacc_out}vacc_paired_phylum_msa_tax_data.csv')
labels, data = _process_column_piechart(vacc_paired_xref, 'class')
class_df = pd.DataFrame(data={'labels':labels, 'counts': data})
class_df.to_csv(f'{vacc_out}vacc_paired_class_msa_tax_data.csv')
labels, data = _process_column_piechart(vacc_paired_xref, 'species')
species_df = pd.DataFrame(data={'labels':labels, 'counts': data})
species_df.to_csv(f'{vacc_out}vacc_paired_species_msa_tax_data.csv')
