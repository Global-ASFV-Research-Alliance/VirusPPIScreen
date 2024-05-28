# These are definitions that allow you to use UniProt's REST API

import os, sys, time, requests, pdb
import pandas as pd
import xml.etree.ElementTree as ET
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar 

def request_rest_data(url, max_retries=4, delay=0.1):
    # input a comlete UniProt REST API url and this will return the data package
    # has built in error codes specific to REST
    retries = 0
    time.sleep(delay)
    while retries < max_retries:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response
            elif response.status_code == 400:
                print(f"Status Code:{response.status_code} Bad request.\nThere is a problem with your input. {url}")
                return None
            elif response.status_code == 404:
                print(f"Status Code:{response.status_code} Not found.]\nThe resource you requested doesn't exist. {url}")
                return None
            elif response.status_code == 410:
                print(f"Status Code:{respose.status_code} Gone.\nThe resource you requested was removed. {url}")
            elif response.status_code == 500:
                print(f"Status Code:{reponse.status_code} Internal server error. Most likely a temporary problem, but if the problem persists please contact us.")
                retries += 1
                time.sleep(delay)
                delay = delay * 2 
                continue
            elif response.status_code == 503:
                print(f"Status code:{reponse.status_code} Service not available. The server is being updated, try again later.")
                retries += 1
                time.sleep(delay)
                delay = delay * 2 
                continue
            else:
                retries += 1
                time.sleep(delay)
                delay = delay * 2 
        except:
            retries += 1
            time.sleep(delay)
            delay = delay * 2 
    print(f"Max retries of {max_retries} reached with url:\n{url}\nExiting...")
    return None

def parse_gene_names(uniprot_txt_data, accession, log_file):
    # this takes the unitprot text data format of a protein and parses out gene names and synonyms
    # currently this colates 'GN'  names from: Name=, Synonyms=, and ORFNames= 
    lines = uniprot_txt_data.split('\n')
    names = []
    f = open(log_file, 'a')
    for line in lines:
        if line.startswith('GN'):
            f.write(f"Original line: {line}\n")
            chunks = line[2:].split(';')
            for chunk in chunks:
                data = chunk.strip().split('=')
                if data[0] == 'Name':
                    try:
                        names.extend(data[1].split(', '))
                    except:
                        breakpoint()
                elif data[0] == 'Synonyms':
                    try:
                        names.extend(data[1].split(', '))
                    except:
                        breakpoint()
                elif data[0] == 'ORFNames':
                    try:
                        names.extend(data[1].split(', '))
                    except:
                        breakpoint()
                else:
                    continue
    name_str = ""
    for name in names:
        name_str += f"{name}, "
    name_str = name_str[:-2]
    f.write(f"parsed names: {name_str}")
    f.write('\n\n')
    f.close()
    return names, name_str

def gene_names_from_accession(prot_accessions=[], log_file="/home/jacob.fenster/scripts/api/uniprot_gene_name_request.log"):
    # input list of Uniprot protein acession numbers and this will request and parse the 
    # gene name and synonmous gene names from it
    os.system(f"> {log_file}")
    summary = pd.DataFrame(columns=['gene_names'])
    summary.index.name = 'uniprot_accession'
    for accession in prot_accessions:
        url = f"https://www.uniprot.org/uniprot/{accession}.txt"
        response = request_rest_data(url)
        if response is None:
            os.system(f'echo "accession {accession} timed out. Skipping...." >> {log_file}')
        elif response.status_code == 200:
            os.system(f'echo "accession {accession} has status code: {response.status_code}" >> {log_file}')
            names, name_str = parse_gene_names(response.text, accession, log_file)
            summary.loc[accession] = {'gene_names': name_str}
        else:
            os.system(f'echo "accession {accession} has status code: {response.status_code}" >> {log_file}')
    return summary

def _parse_aa_seq(uniprot_txt_data):
    lines = uniprot_txt_data.split('\n')
    for i in range(len(lines)):
        if lines[i].startswith('SQ'):
            i += 1
            seq = ''
            while not lines[i].startswith('//') and lines[i].startswith(' '):
                seq += lines[i].strip().replace(' ', '')
                i += 1
    return seq
                

def aa_seq_from_accession(prot_accessions=[], log_file='/home/jacob.fenster/aa.log'):
    os.system(f"> {log_file}")
    summary = pd.DataFrame(columns=['aa_seq'])
    summary.index.name = 'uniprot_accession'
    for accession in prot_accessions:
        url = f"https://www.uniprot.org/uniprot/{accession}.txt"
        response = request_rest_data(url)
        if response is None:
            os.system(f'echo "accession {accession} timed out. Skipping...." >> {log_file}')
        elif response.status_code == 200:
            os.system(f'echo "accession {accession} has status code: {response.status_code}" >> {log_file}')
            seq = _parse_aa_seq(response.text)
            summary.loc[accession] = {'aa_seq': seq}
        else:
            os.system(f'echo "accession {accession} has status code: {response.status_code}" >> {log_file}')
    return summary
