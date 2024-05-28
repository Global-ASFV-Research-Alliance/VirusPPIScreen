import glob, sys, os, pdb
import pandas as pd
import requests
import xml.etree.ElementTree as ET
import time

def request_data_skip(url, max_retries=20, delay=0.10):
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response
            else:
                retries += 1
                time.sleep(delay)
        except:
            retries += 1
            time.sleep(delay)
    print(f"Max retries of {max_retries} reached with url:\n{url}\nSkipping")
    return None

def check_200_request_error(url, response, max_retries=20):
        retries = 0
        while retries < max_retries:
            try:
                root = ET.fromstring(response.text)
                while root.find('./ERROR') != None and retries < max_retries:
                    response = request_data(url)
                    try:
                        root = ET.fromstring(response.text)
                    except: #no automated logic here yet
                        retries += 1
                if retries < max_retries:
                    return response
                else:
                    print(f"Max retries of {max_retries} reached with url below in check_200_request_error:\n{url}\nSkipping")
                    return None
            except: #no automated logic here yet
                retries += 1
                response = request_data(url)

def elink_single_nuccore_taxonomy(search_Ids, api_key, exp_tag, output_dir, 
                                  begin=0, end=None, savepoint=6000, log='/home/jacob.fenster/elink_nuc_tax.log'):
    # This takes a list of accession search Ids and individually queries NCBI for the taxonomy link
    elink_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    if end is None:
        max_count = len(search_Ids)
    else:
        max_count = end
    os.system(f'echo "starting nuccore tax xref for #accessions:{max_count}. start, end index:{begin}, {end}." > {log}')
    xref = {} # crossref dictionary with keys accession and entries taxID
    timei = time.time()
    progress = int(savepoint/10)
    if progress == 0:
        progress = 1
    for index in range(begin, max_count): 
        if index % progress == 0:
            os.system(f'echo "{int(100*((index-begin)/max_count))}% complete" >> {log}')
        # elink chunk
        link_url = f"{elink_base}?dbfrom=nuccore&db=taxonomy&id={search_Ids[index]}&linkname=nuccore_taxonomy&cmd=neighbor&api_key={api_key}"
        link_request = request_data_skip(link_url, max_retries=5)
        if link_request is None:
            os.system(f'echo "ERROR timeout: {search_Ids[index]} accession elink failed due to timeout. Skipping..." >> {log}')
            continue
        link_request = check_200_request_error(link_url, link_request, max_retries=5)
        if link_request is None:
            os.system(f'echo "ERROR timeout: {search_Ids[index]} accession elink failed due to timeout. Skipping..." >> {log}')
            continue
        link_Ids_list = []
        try: # extract linked TaxIDs
            root_link = ET.fromstring(link_request.text)
            for link in root_link.findall(".//Link"):
                link_Ids_list.append(f"{link.find('./Id').text}")
        except:
            os.system(f'echo "ERROR other: {search_Ids[index]} has other error with extracting TaxID. Skipping..." >> {log}')
        if len(link_Ids_list) == 1:
            xref[search_Ids[index]] = link_Ids_list[0]
        elif len(link_Ids_list) > 1:
            os.system(f'echo "ERROR multiple xrefs: {search_Ids[index]}.  Skipping..." >> {log}')
            continue
        else:
            breakpoint()
            os.system(f'echo "ERROR zero xrefs: {len(link_Ids_list)} xrefs for {search_Ids[index]}" >> {log}')
        if index % savepoint  == 0:
            timef = time.time()
            os.system(f'echo "Linked {int(100*((index-begin)/max_count))}% of search query. {timef-timei:0f} seconds elapsed for {savepoint} queries\nSaving progress to {output_dir}{exp_tag}_taxID_xref.csv" >> {log}')
            with open(f"{output_dir}{exp_tag}_taxId_xref.csv", 'w') as f:
                f.write(f"accession,taxID,next_index:{index+1}\n")
                for key in xref:
                    f.write(f"{key},{xref[key]}\n")
    os.system(f'echo "elink_nuccore_taxonomy complete. saving output to {output_dir}{exp_tag}_taxID_xref_final.csv" >> {log}')
    with open(f"{output_dir}{exp_tag}_taxId_xref_final.csv", 'w') as f:
        f.write(f"accession,taxID\n")
        for key in xref:
            f.write(f"{key},{xref[key]}\n")
    return xref 

# data input
api_key = "46f09d13b8f485ab42dabf140b3a09dd9c09"
api_key_usda = '332aca54a1bc5882797d51672572237b6c08'
db = 'nuccore'
output_dir = '/home/jacob.fenster/scripts/api/data/'
acessions_df = pd.read_csv("/home/jacob.fenster/scripts/api/data/20240411_all_virus_searchIds.csv")
acessions = acessions_df['acession'].tolist()
exp_tag= 'troubleshoot_0xref'
log = f'{output_dir}{exp_tag}.log'

def find_error(xml_str):
    root = ET.fromstring(xml_str)
    breakpoint()
    for error in root.find('/ERROR'):
        print(error.text)
elink_single_nuccore_taxonomy(acessions, api_key, exp_tag, output_dir,
                                  begin=0, end=None, savepoint=6000, log=log)

error_xml = '''<?xml version=1.0 encoding=UTF-8 ?>
<!DOCTYPE eLinkResult PUBLIC -//NLM//DTD elink 20101123//EN https://eutils.ncbi.nlm.nih.gov/eutils/dtd/20101123/elink.dtd>
<eLinkResult>
<LinkSet>
	<DbFrom>nuccore</DbFrom>
	<IdList>
	</IdList>
	<ERROR>BLOB ID IS NOT IMPLEMENTED</ERROR>
	<LinkSetDb>
		<DbTo>taxonomy</DbTo>
		<LinkName>nuccore_taxonomy</LinkName>
		<ERROR>NCBI C++ Exception: Error: UNK_MODULE(CException::eInvalid) UNK_FILE, line 18446744073709551615: UNK_FUNC --- </ERROR>
	</LinkSetDb>
</LinkSet>
</eLinkResult>'''
breakpoint()
