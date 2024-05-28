# This script is to filter the MINT database for asfv-asfv ppis

import os, sys, pdb, glob
import pandas as pd
import common.MINT_tools as mint
import common.rbh_tools as rbh
import common.REST_api as rest
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar


asfv_out = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/asfv-asfv_all_homo-heterodimer/"
mint_db_file = "/lustrefs/fadru/projects/asfv-ppi/data/MINT/MINT_db.tsv"

# filter the MINT db by PPIs that are asfv with anything else. asym
asfv_asfv_ppi = mint.filter_mint_sym(mint_db_file, species=['asfv', 'ASFV', 'african swine fever virus', 'African swine fever virus', 
                                         'African swine fever virus ASF', 'African swine fever virus ASFV',
                                         'African swine fever virus, ASFV'],
                                     taxID='10497')
asfv_asfv_ppi.to_csv(f"{asfv_out}MINT_asfv_asfv_ppi.csv")

asfv_asym_ppi = mint.filter_mint_asym(mint_db_file, species=['asfv', 'ASFV', 'african swine fever virus', 'African swine fever virus', 
                                         'African swine fever virus ASF', 'African swine fever virus ASFV',
                                         'African swine fever virus, ASFV'],
                                     taxID='10497')
asfv_asym_ppi.to_csv(f"{asfv_out}MINT_asfv_all_ppi.csv")
