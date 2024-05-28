# This script is to filter the MINT database for vaccinia-vaccinia PPIs
# and crossreference them to the Vaccinia-WR genome via direct UniProt accession number
# or by pulling the amino acid sequence via UniProt REST API and doing a blastp against Vaccinia-WR 
import os, sys, pdb, glob
import pandas as pd
import common.MINT_tools as mint
import common.rbh_tools as rbh
import common.REST_api as rest
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar


vaccinia_out = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/"
vaccinia_wr_proteome = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
vaccinia_wr_genenames = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/vaccinia_all_homo-heterodimer/MINT_vaccinia_genenames.csv"
mint_db_file = "/lustrefs/fadru/projects/asfv-ppi/data/MINT/MINT_db.tsv"

if False:
# filter the MINT db by PPIs that are vaccinia-vaccinia 
    vacc_vacc_ppi = mint.filter_mint_sym(mint_db_file,
                                    species=['Vaccinia', 'vaccinia', 'Vaccinia virus', 'vaccinia virus'], 
                                    taxID='10254')
#vacc_vacc_ppi.to_csv(f"{vaccinia_out}MINT_vaccinia-vaccinia_ppi.csv")

# cross ref exact matches with vaccinia Vaccinia-WR
    crossref_mint = mint._crossref_uniprot(vacc_vacc_ppi,
                                           vaccinia_wr_genenames,
                                           vaccinia_wr_proteome, 
                                           log_file='vaccinia_gold_standards.log', 
                                           temp='/home/jacob.fenster/scripts/mint_uniprot_rbh/tmp/')
    crossref_mint.to_csv(f"{vaccinia_out}MINT_vaccina_ppi_crossref.csv")
breakpoint()
crossref_mint = pd.read_csv(f"{vaccinia_out}MINT_vaccina_ppi_crossref.csv", index_col=0)
vacc_master = pd.read_csv(f"{vaccinia_out}vacc_master_partial8.csv", index_col=0)
vacc_master = mint.load_mint_to_master(crossref_mint, vacc_master)

breakpoint()
vacc_master.to_csv(f"{vaccinia_out}vacc_master_partial9.csv")
