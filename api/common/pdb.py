# These definitions are to access the PDB and parse its file types
import os, sys
import pandas as pd
from Bio.PDB import MMCIF2Dict
#sys.path.append("/home/jacob.fenster/scripts/alphafold-2.3.2/")
#from alphafold.data import parsers
#from alphafold.data import templates
#from alphafold.data import mmcif_parsing
#from alphafold.data import parsers_JF

# parse release dates from mmcif files and output to txt file for rapid future parsing
#def parse_release_dates():
    #pdb_mmcif_dir = "/lustrefs/public/databases/alphafold/pdb_mmcif/mmcif_files/"
    #mmcif_files = glob.glob(f"{pdb_mmcif_dir}*.cif")
    #with open('/home/jacob.fenster/scripts/af_analysis/data/pdb_release_date.txt', 'w') as f:
        #for file in mmcif_files:
            #mmcif_dict = MMCIF2Dict.MMCIF2Dict(file)
            #release_date = mmcif_parsing.get_release_date(mmcif_dict)
            #pdb_id = mmcif_dict['_entry.id'][0]
            #f.write(f"{pdb_id}:{release_date}\n")

# parse pdb entries.idx file 
def parse_line(line):
    return [field.strip() for field in line.split('\t')]

def entries_idx_parse(input_file):
    # to parse the pdb entries.idx file from the PDB. This ended up having the upload dates, not the release datesc
    # https://files.wwpdb.org/pub/pdb/derived_data/index/ 
    # 
    parsed_data = []
    with open(input_file, 'r') as file:
        next(file)
        next(file)
        for line in file:
            if line.strip():
                parsed_data.append(parse_line(line))

    columns = ['IDCODE', 'HEADER', 'ACCESSION DATE', 'COMPOUND', 'SOURCE', 'AUTHOR LIST', 'RESOLUTION', 'EXPERIMENT TYPE (IF NOT X-RAY)']
    df = pd.DataFrame(parsed_data, columns=columns)
    return df
