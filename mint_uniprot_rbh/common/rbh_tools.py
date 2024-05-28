# This contains functions to conduct reciprocal best hit analysis using BLASTp
import os, sys, subprocess, pdb
import pandas as pd
import numpy as np

# Global paths to compiled NCBI-BLAST software
blastp = "/home/jacob.fenster/scripts/rbh/common/ncbi-blast/bin/blastp"
makeblastdb = "/home/jacob.fenster/scripts/rbh/common/ncbi-blast/bin/makeblastdb"
# Global type dictionary for blastp format specifiers 
outfmt_type = {'qseqid': str, 'sseqid': str, 'evalue': float, 
               'pident': float, 'qcovs': float, 'slen': int, 
               'qlen': int, 'length': int, 'qstart': int, 
               'qend': int, 'sstart': int, 'send': int,
               'bitscore': float, 'qseq': str, 'sseq': str}

def _makeblastdb(fasta_file, output_dir, dbtype='prot',  make_blastdb=makeblastdb):
    """ 
    Make a blastdb from an input fasta file. Returns the path to the database or None if failed 
    does not attempt to parse the fasta descriptions 
    """
    basename = os.path.basename(fasta_file)
    if os.path.exists(f"{output_dir}{basename}.pjs"):
        return f"{output_dir}{basename}"
    else:
        command = f'{makeblastdb} -in "{fasta_file}" -blastdb_version 5 -out {output_dir}{basename} -dbtype {dbtype}' 
        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        if result.returncode == 0:
            return f"{output_dir}{basename}" 
        else:
            print(f"ERROR: {result.stderr}")
            print(f"ARGS: {result.args}")
            print(f"STDOUT: {result.stdout}")
            return None

def _blastp_search(query_fasta, database, output_dir='-', 
                   outfmt='10 qseqid sseqid evalue pident qcovs slen qlen length qstart qend sstart send',
                   sorthits=3, evalue=.01, blastp=blastp,
                   max_target_seqs=100,
                   num_threads=1, mt_mode=1):
    """
    Conducts a blastp search of the query fasta against the input database 
    if no output_dir specified it outputs the results to stdout 
    """
    if output_dir == '-':
        basename = ''
    else:
        basename = f"{os.path.splitext(os.path.basename(query_fasta))[0]}__{os.path.basename(database)}.tsv"
    command = f'{blastp} -query "{query_fasta}" -db "{database}" -out "{output_dir}{basename}" -evalue {evalue} \
-max_target_seqs {max_target_seqs} -sorthits {sorthits} -outfmt "{outfmt} -mt_mode {mt_mode} -num_threads {num_threads}"' 
    result = subprocess.run(command, shell=True, text=True, capture_output=True)
    if result.returncode == 0:
        if output_dir == '-':
            return result.stdout, outfmt
        else:
            return f"{output_dir}{basename}", outfmt 
    else:
        print(f"ERROR {result.returncode}: {result.stderr}")
        print(f"ARGS: {result.args}")
        return None, None

def _blastp_to_df(query_fasta, subject_fasta, db_out, 
           outfmt='10 qseqid sseqid evalue pident qcovs slen qlen length qstart qend sstart send',
           evalue=0.01):
    """
    does a blastp search of the query versus subject and loads the output into a dataframe
    """
    # generate df with appropriate columns 
    columns = []
    # load column titles
    for entry in outfmt.split():
        if entry == '10':
            continue
        columns.append(entry)
    # generate database and search 
    subject_db = _makeblastdb(subject_fasta, db_out)
    blast_out, outfmt = _blastp_search(query_fasta, subject_db, outfmt=outfmt, evalue=evalue)
    # see if no hits or if incorrect outfmt
    if blast_out == '':
        print(f"No hits returned searching subject {subject_fasta} with query {query_fasta} with evalue threshold of {evalue}")
        return pd.DataFrame(columns=columns)
    elif len(blast_out.split('\n')[0].split(',')) != len(columns):
        print(f"ERROR with outfmt: {outfmt}")
        return None
    # build dataframes with blastp output
    elif blast_out == None:
        return None
    else:
        rows = []
        for line in blast_out.split('\n'):
            entries = line.split(',')
            row = {}
            if entries != ['']: # artifact with last line
                for title, entry in zip(columns, entries):
                    row[title] = entry
                rows.append(row)
        blast_df = pd.DataFrame(data=rows, columns=columns)
        # deal with type formatting
        for column in columns:
            try:
                blast_df[column] = blast_df[column].astype(outfmt_type[column])
            except KeyError:
                print(f"KeyError: column {column} not in outfmt_type dictionary")
        return blast_df

def _hit_overlap(qstarts, qends):
    """ 
    calculates the overlap length between current hit and previous hits 
    qstarts and qends are list of indicies of the query start and query end index
    """
    qstart = max(qstarts)
    qend = min(qends)
    if qstart > qend:
        return 0
    else:
        return qend - qstart 

def _best_hit_procedure(blast_df,
                        hit_len_tr=50, qcovs_tr=50, hcovs_tr=50):
    """
           max_target_seqs=100,
    This is a forward/reverse best hit algorithm mirrored from Ian et al 2021. 
    blast_df is the blast output in a dataframe that has gone through evalue prefiltering
    and has hit coverage, hcovs, precomputed 
    returns a list of best hit dataframes for each unique query
    """
    # calculate hit coverage
    blast_df['hcovs'] = ((blast_df['send'] - blast_df['sstart']) / blast_df['slen']) * 100
    # filter blast hit 
    blast_df = blast_df[(blast_df['length'] > hit_len_tr) | (blast_df['qcovs'] > qcovs_tr)]
    blast_df = blast_df[(blast_df['hcovs'] > hcovs_tr) | (blast_df['qcovs'] > qcovs_tr)]
    # 
    qseqid_set  = set(blast_df['qseqid'].tolist())
    best_hits = [] #list of best hit dataframes
    for qseqid in qseqid_set:
        seqid_df = blast_df[blast_df['qseqid'] == qseqid]
        seqid_pident = seqid_df.sort_values(by='pident', ascending=False)
        seqid_evalue = seqid_df.sort_values(by='evalue', ascending=True)
        # process rank by pident
        first = True
        best_id = max(seqid_pident['pident'].tolist())
        hits_pident = pd.DataFrame(columns=seqid_df.columns)
        for index in seqid_pident.index:
            if first:
                hits_pident.loc[index] = seqid_pident.loc[index]
                first = False
            else:
                qstarts = [seqid_pident.loc[index, 'qstart']] + hits_pident['qstart'].tolist()
                qends = [seqid_pident.loc[index, 'qend']] + hits_pident['qend'].tolist()
                overlap = _hit_overlap(qstarts, qends)
                if overlap > 0.25 * seqid_pident.loc[index, 'length'] and seqid_pident.loc[index, 'pident'] > best_id - 10:
                    hits_pident.loc[index] = seqid_pident.loc[index]
                elif overlap < 0.25 * seqid_pident.loc[index, 'length'] and seqid_pident.loc[index, 'pident'] > best_id - 20:
                    hits_pident.loc[index] = seqid_pident.loc[index]
        # process rank by evalue
        first = True
        best_evalue = min(seqid_evalue['evalue'].tolist())
        hits_evalue = pd.DataFrame(columns=seqid_df.columns)
        for index in seqid_evalue.index:
            if first:
                hits_evalue.loc[index] = seqid_evalue.loc[index]
                first = False
            else:
                qstarts = [seqid_evalue.loc[index, 'qstart']] + hits_evalue['qstart'].tolist()
                qends = [seqid_evalue.loc[index, 'qend']] + hits_evalue['qend'].tolist()
                overlap = _hit_overlap(qstarts, qends)
                if overlap > 0.25 * seqid_evalue.loc[index, 'length'] and seqid_evalue.loc[index, 'evalue'] < best_evalue * 1e5:
                    hits_evalue.loc[index] = seqid_evalue.loc[index]
                elif overlap < 0.25 * seqid_evalue.loc[index, 'length']:
                    hits_evalue.loc[index] = seqid_evalue.loc[index]
        common_indicies = hits_evalue.index.intersection(hits_pident.index)
        intersection = hits_evalue.loc[common_indicies]
        # can rename the index to the query hit id? 
        best_hits.append(intersection)
    return best_hits

def rbh(proteomeA, proteomeB, db_out, 
       hit_len_tr=50, qcovs_tr=50, hcovs_tr=50,   
       outfmt='10 qseqid sseqid evalue pident bitscore qcovs slen qlen length qstart qend sstart send',
       evalue=0.01):
    """ 
    Conducts Ian humphrey's et al. 2021 reciprocal best hit orthologue finder
    hit_len_tr is the base length that the length of alignment must be to pass filter
    qcovs_tr is the query coverage threshold, percentage, that must be met to pass filter
    hcovs_tr is the hit coverage threshold, percentage, that must be met to pass filter
    """
    # conduct the blastp search

    a_query_b = _blastp_to_df(query_fasta=proteomeA, subject_fasta=proteomeB, db_out=db_out, outfmt=outfmt, evalue=evalue)
    b_query_a = _blastp_to_df(query_fasta=proteomeB, subject_fasta=proteomeA, db_out=db_out, outfmt=outfmt, evalue=evalue)
    # calculate hit coverage, hcovs
    a_query_b_hits = _best_hit_procedure(a_query_b)
    b_query_a_hits = _best_hit_procedure(b_query_a)
    # There is room to optimize the loop here with dictionaries or listso
    rbh = pd.DataFrame(columns=['query', 'subject', 'pair'])
    for hit_f in a_query_b_hits:
        query, subject_ids = hit_f['qseqid'].tolist()[0], set(hit_f['sseqid'].tolist())
        for hit_r in b_query_a_hits:
            subject, query_ids = hit_r['qseqid'].tolist()[0], set(hit_r['sseqid'].tolist())
            if subject in subject_ids and query in query_ids:  
                new_row = pd.DataFrame({'query': [query], 'subject': [subject], 'pair': [f"{query}:{subject}"]})
                rbh = pd.concat([rbh, new_row], ignore_index=True)

    return rbh

def fbh(queries, proteome, db_out, 
       hit_len_tr=50, qcovs_tr=50, hcovs_tr=50,   
       outfmt='10 qseqid sseqid evalue pident bitscore qcovs slen qlen length qstart qend sstart send',
       evalue=0.01):
    """
    this does a forward best hit of the queries (fasta format) against the proteome (fasta format)
    values and thresholds taken from Humphreys et al 2021. 
    """
    blast_df = _blastp_to_df(query_fasta=queries, subject_fasta=proteome, db_out=db_out, outfmt=outfmt, evalue=evalue)
    hits = _best_hit_procedure(blast_df)
    return hits

