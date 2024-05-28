import os, sys, glob, subprocess, time, datetime, pdb
import pandas as pd
import numpy as np
import analysis.AFppi as ppi
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar
import pmsa_tools as pmsa
sys.path.append("/home/jacob.fenster/scripts/alphafold-2.3.2/")
from alphafold.data import parsers
from alphafold.data import templates
from alphafold.data import mmcif_parsing
from alphafold.data import parsers_JF
from Bio.PDB import MMCIF2Dict

def hhsearch_to_TemplateHit(
        protein='',
        af_output_dir='',
        tmp='/home/jacob.fenster/tmp/', 
        pdb70="/lustrefs/public/databases/alphafold/pdb70/", 
        hhsearch='hhsearch', 
        reformat= '/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl',
    ):
    """
    This converts a uniref90 MSA sto file to a3m
    does a hhsearch of the pdb70 db, parses it with alphafold.parsers to convert into a list of TemplateHit objects 
    returns the query sequence of the input uniref90 sto file and the list of TemplateHits
    """
    os.makedirs(tmp, exist_ok=True)
    msas =  pmsa.read_AF_msas(protein, af_output_dir)
    if msas is not None:
        for msa in msas:
            if os.path.basename(msa) == "uniref90_hits.sto":
                uniref90_sto = msa
    else:
        breakpoint()
    reformat_command = f"{reformat} sto a3m -v 0 {uniref90_sto} {tmp}{protein}-uniref90_hits.a3m"
    hhsearch_command = f"{hhsearch} -i {tmp}{protein}-uniref90_hits.a3m -d {pdb70}pdb70 -o {tmp}{protein}-pdb70hits.hhr -maxres 1000000"
    try:
        reformat_result = subprocess.run(reformat_command, shell=True, capture_output=True)
        descs, seqs = bpar.read_fasta(f"{tmp}{protein}-uniref90_hits.a3m")
        query_sequence = seqs[0].replace('-', '').replace('.', '')
        hhsearh_result = subprocess.run(hhsearch_command, shell=True, capture_output=True)
        with open(f"{tmp}{protein}-pdb70hits.hhr", 'r') as f:
            hrr_string = f.read()
    except:
        breakpoint()
    # this returns a list of TemplateHit objects
    template_hits = parsers.parse_hhr(hrr_string)
    return query_sequence, template_hits

def parse_file_hrr(file_path):
    with open(file_path, 'r') as f:
        hrr_string = f.read()
    template_hits = parsers.parse_hhr(hrr_string)
    return template_hits

def parse_file_hrr_JF(file_path):
    with open(file_path, 'r') as f:
        hrr_string = f.read()
    template_hits = parsers_JF.parse_hhr_JF(hrr_string)
    return template_hits

def assess_hhsearch_hit_jf(
    hit: parsers.TemplateHit,
    hit_pdb_code: str,
    query_sequence: str,
    release_dates,
    release_date_cutoff: datetime.datetime,
    id_tr: float=95.0,
    max_subsequence_ratio: float=.95,
    min_subsequence_ratio: float=0.1):
  """Determines if template is valid (without parsing the template mmcif file).

  Args:
    hit: HhrHit for the template.
    hit_pdb_code: The 4 letter pdb code of the template hit. This might be
      different from the value in the actual hit since the original pdb might
      have become obsolete.
    query_sequence: Amino acid sequence of the query.
    release_dates: Dictionary mapping pdb codes to their structure release
      dates.
    release_date_cutoff: Max release date for desired AF params
    max_subsequence_ratio: Exclude any exact matches with this much overlap.
    min_subsequence_ratio: minimum ratio of identical sequences allowed to be called a partial_duplicate
    id_tr: identities threshold to call a duplicate hit
  Returns:
      assessment dictionary which has structure {'included_af_params': Bool, 'duplicate': Bool, 'partial_duplicate': Bool}
###################################################################
    JF edits: trying to change this function so that it returns:
    - if template is included in the run model parameters
    - if there is an exact template match and the ID of the exact match
    - Future - tell us if there is also an exact match with pdb heteromer structure
  """
  aligned_cols = hit.aligned_cols
  identities = hit.identities
  #align_ratio = aligned_cols / len(query_sequence)

  template_sequence = hit.hit_sequence.replace('-', '')
  template_len = len(template_sequence)
  length_ratio = float(len(template_sequence)) / len(query_sequence)

  # Check whether the template is a large subsequence or duplicate of original
  # query. This can happen due to duplicate entries in the PDB database.
  duplicate = (identities > id_tr and length_ratio > max_subsequence_ratio)
  partial_duplicate = (identities > id_tr and length_ratio > min_subsequence_ratio)
  #need to add some % ID similarity to find strain variants. use 95% ID cuttoff 
  included_af_params = (not templates._is_after_cutoff(hit_pdb_code, release_dates, release_date_cutoff))

  assessment = {} #to hold hit assessment results
  if duplicate:
    assessment['duplicate'] = True
    assessment['partial_duplicate'] = True
    if included_af_params:
        assessment['included_af_params'] = True
    else:
        assessment['included_af_params'] = False 
  elif partial_duplicate:
    assessment['duplicate'] = False
    assessment['partial_duplicate'] = True
    if included_af_params:
        assessment['included_af_params'] = True
    else:
        assessment['included_af_params'] = False 
  else:
      #potential here to include coverage and %ID to template
    assessment['duplicate'] = False
    assessment['partial_duplicate'] = False
    assessment['included_af_params'] = False
  assessment['identities'], assessment['template_length'] = identities, template_len
  return assessment

# parse release dates from mmcif files and output to txt file for rapid future parsing
def parse_release_dates():
    pdb_mmcif_dir = "/lustrefs/public/databases/alphafold/pdb_mmcif/mmcif_files/"
    mmcif_files = glob.glob(f"{pdb_mmcif_dir}*.cif")
    with open('/home/jacob.fenster/scripts/af_analysis/data/pdb_release_date.txt', 'w') as f:
        for file in mmcif_files:
            mmcif_dict = MMCIF2Dict.MMCIF2Dict(file)
            release_date = mmcif_parsing.get_release_date(mmcif_dict)
            pdb_id = mmcif_dict['_entry.id'][0]
            f.write(f"{pdb_id}:{release_date}\n")

# parse pdb entries.idx file 
def parse_line(line):
    return [field.strip() for field in line.split('\t')]

def idx_parse(input_file):
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

# this script is attempting to determine if PDB templates existed as input for the original AF model or the feature input
# starting by looking at 12L__A40R. 12L is D12L and is 287 aa long 
output_dir = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/af_template_analysis/"
tmp = f"{output_dir}tmp/" 
afmult_vaccinia_out = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/AFmult/"
release_dates = templates._parse_release_dates("/lustrefs/fadru/projects/asfv-ppi/data/pdb_release_dates.txt")
max_template_date = "2018-04-30" #AlphaFold Multimer parameters training
date_format = "%Y-%m-%d"
release_date_cutoff = datetime.datetime.strptime(max_template_date, date_format)

template_analysis = pd.DataFrame(columns=['duplicate', 'partial_duplicate','included_af_params',
                                          'query_length', 'max_id', 'max_id_length', 'max_length', 'max_length_id'])
mccraith_proteins = ["D12L", "A40R"]

mccraith_afmult_fastas = glob.glob(f"/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/*.fasta")
mccraith_proteins = set()
pairs = [] # holds McCraith PPI pairs parsed from AF mult run
for fasta in mccraith_afmult_fastas: 
    descs, seqs = bpar.read_fasta(fasta)
    pairs.append((descs[0].split('__')[0], descs[1].split('__')[0]))
    mccraith_proteins.add(descs[0].split('__')[0]), mccraith_proteins.add(descs[1].split('__')[0])

# load pre computed hsearch .hrr files 
for protein in mccraith_proteins:
    template_hits = parse_file_hrr_JF(f"{tmp}{protein}-pdb70hits.hhr")
    descs, seqs = bpar.read_fasta(f"{tmp}{protein}-uniref90_hits.a3m")
    query_sequence = seqs[0].replace('-', '').replace('.', '').upper()
    query_len = len(query_sequence)
    af_params, duplicate, partial_duplicate = False, False, False
    identities_list, template_len_list = [], []
    hit_identities_list, hit_template_len_list = [], []
    #cycle through all the hits and record if any are duplicate or partial duplicates. Record metrics from those hits
    for template_hit in template_hits:
        id, chain = templates._get_pdb_id_and_chain(template_hit)
        assessment = assess_hhsearch_hit_jf(template_hit, id, query_sequence, release_dates, release_date_cutoff, id_tr=95.0, 
                                            max_subsequence_ratio=0.95, min_subsequence_ratio=0.1)
        identities_list.append(assessment['identities']), template_len_list.append(assessment['template_length'])
        if assessment['included_af_params']:
            af_params = True
        if assessment['duplicate']:
            duplicate = True
        if assessment['partial_duplicate']:
            partial_duplicate = True
        if assessment['duplicate'] or assessment['partial_duplicate']:
            hit_identities_list.append(assessment['identities']), hit_template_len_list.append(assessment['template_length'])
    if assessment['duplicate'] or assessment['partial_duplicate']:
        max_id = max(hit_identities_list)
        max_id_len = hit_template_len_list[hit_identities_list.index(max_id)]
        max_len = max(hit_template_len_list)
        max_len_id = hit_identities_list[hit_template_len_list.index(max_len)]
    else:
        max_id = max(identities_list)
        max_id_len = template_len_list[identities_list.index(max_id)]
        max_len = max(template_len_list)
        max_len_id = identities_list[template_len_list.index(max_len)]
    template_analysis.loc[protein] = {'duplicate': duplicate, 'partial_duplicate': partial_duplicate, 'included_af_params': af_params,
                                      'query_length': query_len, 'max_id': max_id, 'max_id_length': max_id_len, 'max_length': max_len, 
                                      'max_length_id': max_len_id}

template_analysis.to_csv(f"{output_dir}McCraith2000_individual_template.csv")
breakpoint()





if False: #to compute the hhsearch
    for protein in mccraith_proteins:
        for template_hit in template_hits:
            id, chain = templates._get_pdb_id_and_chain(template_hit)
            assessment = assess_hhsearch_hit_jf(template_hit, id, query_sequence, release_dates, release_date_cutoff, id_tr=95.0,
                                                max_subsequence_ratio=0.95, min_subsequence_ratio=0.1)
            if assessment['included_af_params']:
                af_params = True
            if assessment['duplicate']:
                duplicate = True
            if assessment['partial_duplicate']:
                partial_duplicate = True
        template_analysis.loc[protein] = {'included_af_params': af_params, 'duplicate': duplicate, 'partial_duplicate': partial_duplicate}
    template_analysis.to_csv(f"{output_dir}McCraith2000_individual_template.csv")
        query_sequence, template_hits = hhsearch_to_TemplateHit(protein, af_output_dir=afmult_vaccinia_out, tmp=tmp)
        af_params, duplicate, partial_duplicate = False, False, False
#pa = pd.read_csv(f"{output_dir}McCraith2000_individual_template.csv", index_col=0)
#pair_analysis = pd.DataFrame(columns=['protA', 'protA_duplicate', 'protA_partial_duplicate', 'protA_included', 
                                      'protB', 'protB_duplicate', 'protB_partial_duplicate', 'protB_included',
                                      'both_duplicate', 'one_duplicate', 'mixed_duplicate_partial', 'partial_duplicate', 'None',
                                      'both_included', 'one_included', 'none_included'])
if False:
    for pair in pairs:
        pair_name = f"{pair[0]}__{pair[1]}"
        protA, protB = pair[0], pair[1]
        pair_analysis.loc[pair_name] = {'protA': protA, 'protA_duplicate': pa.loc[protA, 'duplicate'], 'protA_partial_duplicate': pa.loc[protA, 'partial_duplicate'], 
                                        'protA_included': pa.loc[protA, 'included_af_params'],
                                        'protB': protB, 'protB_duplicate': pa.loc[protB, 'duplicate'], 'protB_partial_duplicate': pa.loc[protB, 'partial_duplicate'], 
                                        'protB_included': pa.loc[protB, 'included_af_params'],
                                        'both_duplicate': (pa.loc[protA, 'duplicate'] and pa.loc[protB, 'duplicate']), 
                                        'one_duplicate':(pa.loc[protA, 'duplicate'] or  pa.loc[protB, 'duplicate']), 
                                        'mixed_duplicate_partial': (pa.loc[protA, 'duplicate'] or  pa.loc[protB, 'duplicate']) and (pa.loc[protA, 'partial_duplicate'] or  pa.loc[protB, 'partial_duplicate']), 
                                        'partial_duplicate': (pa.loc[protA, 'partial_duplicate'] or  pa.loc[protB, 'partial_duplicate']),
                                        'None': (not pa.loc[protA, 'duplicate'] and not pa.loc[protB, 'duplicate'] and not pa.loc[protA, 'partial_duplicate'] and not pa.loc[protB, 'partial_duplicate']),
                                        'both_included': (pa.loc[protA, 'included_af_params'] and pa.loc[protB, 'included_af_params']),
                                        'one_included': (pa.loc[protA, 'included_af_params'] or pa.loc[protB, 'included_af_params']),
                                        'none_included': (not pa.loc[protA, 'included_af_params'] and not pa.loc[protB, 'included_af_params'])}
    pair_analysis.to_csv(f"{output_dir}McCraith2000_paired_template.csv")


