# This script is to cross reference each viral genome in the NCBI database with its taxID and other metadata
#
import os, sys, glob, pdb, time
import pandas as pd
import numpy as np
import common.Entrez_NCBI as entrez

# data input
api_key = "46f09d13b8f485ab42dabf140b3a09dd9c09"
api_key_usda = '332aca54a1bc5882797d51672572237b6c08'
db = 'nuccore'
output_dir = '/home/jacob.fenster/scripts/api/data/'
acessions_df = pd.read_csv("/home/jacob.fenster/scripts/api/data/20240411_all_virus_searchIds.csv")
acessions = acessions_df['acession'].tolist()

tranche, size, intervals, api_key, exp_tag = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
begin = (tranche-1)*int(size/intervals)
if tranche == intervals:
    end = size
else:
    end = (tranche)*int(size/intervals)
xref = entrez.elink_single_nuccore_taxonomy(acessions, api_key, exp_tag, output_dir, begin=begin, end=end, savepoint=6000, log=f'{output_dir}{exp_tag}_elink.log')
