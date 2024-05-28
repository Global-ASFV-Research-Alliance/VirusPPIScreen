import requests, time, sys, ast, os, pdb
import xml.etree.ElementTree as ET
import pandas as pd

def request_data(url, max_retries=100, delay=0.10):
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
    print(f"Max retries of {max_retries} reached with url:\n{url}\nExiting")
    sys.exit()


def request_data_skip(url, max_retries=20, delay=0.10):
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url)
            if response.status_code == 200:
                return response
            else:
                time.sleep(delay*(2**retries))
                retries += 1
        except:
            time.sleep(delay*(2**retries))
            retries += 1
    print(f"Max retries of {max_retries} reached with url:\n{url}\nSkipping")
    return None

def check_200_request_error(url, response, max_retries=20):
    """
    This one sees if there are any elements named "ERROR" when the request status is 200. 
    Then does a exponential delayed request 
    """
    retries = 0
    while retries < max_retries:
        try:
            root = ET.fromstring(response.text)
            while len(root.findall('.//ERROR')) != 0 and retries < max_retries:
                response = request_data_skip(url)
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
            response = request_data_skip(url)

def get_element_text(element, tag, nullreturn):
    """Safely get the text from an XML element."""
    sub_element = element.find(tag)
    if sub_element is not None:
        return sub_element.text
    return nullreturn

def parse_nuccore_INSDSeq(entry):
    main = {} #dictionary to hold header information
    source = {} #dictionary to hold sequence feature 'source' information
    main['acession'] = get_element_text(entry, "INSDSeq_primary-accession", None)
    main['locus'] = get_element_text(entry, "INSDSeq_locus", None)
    main['moltype'] = get_element_text(entry, "INSDSeq_moltype", None)
    main['seqlength'] = get_element_text(entry, "INSDSeq_length", None)
    main['topology'] = get_element_text(entry, "INSDSeq_topology", None)
    main['strandedness'] = get_element_text(entry, "INSDSeq_strandedness", None)
    main['definition'] = get_element_text(entry, "INSDSeq_definition", None)
    main['create_date'] = get_element_text(entry, "INSDSeq_create-date", None)
    main['update_date'] = get_element_text(entry, "INSDSeq_create-date", None)
    main['source'] = get_element_text(entry, "INSDSeq_update-date", None)
    main['organsim'] = get_element_text(entry, "INSDSeq_organism", None)
    main['tax_lineage'] = get_element_text(entry, "INSDSeq_taxonomy", None)
    features = entry.findall('.//INSDFeature')
    if features is not None:
        for feature in features:
            if get_element_text(feature, "INSDFeature_key", None) == 'source':
                feature_quals = feature.findall('.//INSDQualifier')
                if feature_quals is not None:
                    for feature_qual in feature_quals:
                        name = get_element_text(feature_qual, "INSDQualifier_name", None)
                        value = get_element_text(feature_qual,"INSDQualifier_value", None)
                        if name in source: #allow ability to append multiple entries
                            source[name] += f"; {value}"
                        else: 
                            source[name] = value
    return main, source

def esearch_id_list(query, db, api_key, exp_tag, output_dir, retmax=10000):
    #this conducts a search and saves the IDs to a output file. 
    #written for large IDs so that I can fetch manually
    print("esearch_id_list started\n")
    base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    url = f"{base}?db={db}&term={query}&idtype=acc&api_key={api_key}&usehistory=y"
    request = request_data(url)
    request = check_200_request_error(url, request)
    try: #try and extract the Count information
        root = ET.fromstring(request.text)
        max_count, QueryKey, WebEnv = int(root.find('./Count').text), root.find('./QueryKey').text, root.find('./WebEnv').text
    except: #this is here to troublshoot. Nothing automated to handle exceptions yet
        pdb.set_trace()
    print(f"There are {str(max_count)} results retured for query:\n {query.replace('+', ' ')}")
    search_Ids = []
    for retstart in range(0, max_count, retmax):
        history_url = f"{base}?db={db}&idtype=acc&query_key={QueryKey}&WebEnv={WebEnv}&retstart={retstart}&retmax={retmax}&usehistory=y&api_key={api_key}"
        history_request = request_data(history_url)
        history_request = check_200_request_error(history_url, history_request)
        try: #extract search acession numbers 
            root_search = ET.fromstring(history_request.text)
            for Id in root_search.findall(".//Id"):
                search_Ids.append(Id.text)
        except:
            pdb.set_trace()
        if retstart % 120000 == 0:
            print(f"Retrieved {int(100*((retstart+retmax)/max_count))}% of search acession numbers")
            print(f"Next restart is {retstart+retmax}\nSaving progress...")
            with open(f"{output_dir}/{exp_tag}_searchIds.csv", 'w') as f:
                f.write(f"acession,Next_retstart:{retstart+retmax},\n")
                for Id in search_Ids:
                    f.write(f'{Id},\n')   
    with open(f"{output_dir}/{exp_tag}_searchIds.csv", 'w') as f:
        f.write(f"acession,Next_retstart:{retstart+retmax},\n")
        for Id in search_Ids:
            f.write(f'{Id},\n')   
    print("esearch_id_list complete\n")
    return search_Ids

def elink_nuccore_taxonomy(search_Ids, tax_master_list, api_key, exp_tag, output_dir, begin=0, retmax=200):
    #this will assemble a df with every taxID and corresponding lineage information 
    print("elink_nuccore_taxonomy started\n")
    elink_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi'
    max_count = len(search_Ids)
    for retstart in range(begin, max_count, retmax): 
        if retstart % 60000 == 0:
            print(f"Elink {int(100*(retstart/max_count))}% complete")
        # elink chunk
        Ids = ''
        try:
            chunk = search_Ids[retstart:retstart+retmax]
            for Id in chunk:
                Ids += f"{Id},"
        except IndexError:
            chunk = search_Ids[retstart:]
            for Id in chunk:
                Ids += f"{Id},"
        Ids = Ids.strip(',')
        link_url = f"{elink_base}?dbfrom=nuccore&db=taxonomy&id={Ids}&linkname=nuccore_taxonomy&cmd=neighbor&api_key={api_key}"
        link_request = request_data(link_url)
        link_request = check_200_request_error(link_url, link_request)
        link_Ids_list = []
        try: # extract linked TaxIDs
            root_link = ET.fromstring(link_request.text)
            for link in root_link.findall(".//Link"):
                link_Ids_list.append(f"{link.find('./Id').text}")
        except:
            pdb.set_trace()
        #filter for new TaxIDs to reduce fetch redundancy
        for taxId in link_Ids_list:
            if taxId not in tax_master_list:
                tax_master_list.append(taxId)
        if retstart % 30000 == 0:
            print(f"Linked {int(100*((retstart+retmax)/max_count))}% of search query. Next retstart: {retstart+retmax}\nSaving progress...")
            with open(f"{output_dir}/{exp_tag}_taxId_list.csv", 'w') as f:
                f.write(f"TaxID,Next_retstart:{retstart+retmax}\n")
                for taxId in tax_master_list:
                    f.write(f"{taxId},\n")
    print("elink_nuccore_taxonomy complete\n")
    return tax_master_list
    
def elink_single_nuccore_taxonomy(search_Ids, api_key, exp_tag, output_dir, 
                                  begin=0, end=None, savepoint=6000, log='/home/jacob.fenster/elink_nuc_tax.log'):
    """
    search_Ids: list []
    """
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
            os.system(f'echo "{100*((index-begin)/max_count):.1f}% complete" >> {log}')
        # elink chunk
        link_url = f"{elink_base}?dbfrom=nuccore&db=taxonomy&id={search_Ids[index]}&linkname=nuccore_taxonomy&cmd=neighbor&api_key={api_key}"
        breakpoint()
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
            os.system(f'echo "ERROR other: {search_Ids[index]} has other error with extracting TaxID. Skipping. xml below" >> {log}')
            os.system(f'echo "{link_request.text}" >> {log}')
        if len(link_Ids_list) == 1:
            xref[search_Ids[index]] = link_Ids_list[0]
        elif len(link_Ids_list) > 1:
            os.system(f'echo "ERROR multiple xrefs: {search_Ids[index]}. xml below" >> {log}')
            os.system(f'echo "{link_request.text}" >> {log}')
            link_string = ''
            for link in link_Ids_list:
                link_string += f'{link};'
            link_string = link_string[:-1]
            xref[search_Ids[index]] = link_string
        else:
            os.system(f'echo "ERROR zero xrefs: {len(link_Ids_list)} xrefs for {search_Ids[index]}. xml below" >> {log}')
            os.system(f'echo "{link_request.text}" >> {log}')
        if index % savepoint  == 0:
            timef = time.time()
            os.system(f'echo "Linked {100*((index-begin)/max_count):.1f}% of search query. {timef-timei:0f} seconds elapsed for {savepoint} queries\nSaving progress to {output_dir}{exp_tag}_taxID_xref.csv" >> {log}')
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

def efetch_taxonomy(tax_master_list, tax_df, api_key, output_dir, exp_tag, begin=0, retmax=100):
    print("efetch_taxonomy started\n")
    efetch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    # efetch taxonomy xml files
    max_count = len(tax_master_list)
    rettype, retmode = 'native', 'xml'
    for retstart in range(begin, max_count, retmax):
        taxId_fetch = ''
        try:
            taxId_slice = tax_master_list[retstart:retstart+retmax]
        except IndexError:
            taxId_slice = tax_master_list[retstart:]
        for taxId in taxId_slice: #assemble API string
            taxId_fetch += f"{taxId},"
        taxId_fetch = taxId_fetch.strip(',')
        fetch_url = f"{efetch_base}?db=taxonomy&id={taxId_fetch}&rettype={rettype}&rettmode={retmode}&api_key={api_key}"
        fetch_request = request_data(fetch_url)
        fetch_request = check_200_request_error(fetch_url, fetch_request)
        try: # extract taxonomy information
            root_tax = ET.fromstring(fetch_request.text)
            for entry in root_tax.findall("Taxon"): #step through taxonomy entries
                taxId = entry.find("TaxId").text
                entry_rank = entry.find("Rank").text 
                ParentTaxId = entry.find("ParentTaxId").text
                tax_df.loc[taxId, 'entry_rank'] = entry_rank
                tax_df.loc[taxId, 'ParentTaxId'] = ParentTaxId
                tax_df.loc[taxId, entry_rank] = entry.find("ScientificName").text
                for lineage in entry.findall("LineageEx"):
                    for taxon in lineage.findall("Taxon"):
                        rank = taxon.find("Rank").text
                        tax_df.loc[taxId, rank] = taxon.find("ScientificName").text
        except:
            pdb.set_trace()
        if retstart % 1000 == 0:
            print(f"Taxonomy fetch {int(100*(retstart/max_count))}% complete\nSaving Progress...")
            with open(f'{output_dir}/{exp_tag}_taxonomy_nextretstart.txt', 'w') as f:
                f.write(f'Next_taxonomy_retstart:{retstart+retmax}')
            tax_df.to_csv(f"{output_dir}/{exp_tag}_taxonomy.csv")
    print("efetch_taxonomy complete\n")
    return tax_df

def efetch_metatadata_nuccore(nuccore_meta_df, search_Ids, api_key, output_dir, exp_tag, begin=0, retmax=100):
    print("efetch_metatadata_nuccore started\n")
    efetch_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
    # fetch metadata
    max_count = len(search_Ids)
    for retstart in range(begin, max_count, retmax):
        Id_fetch = ''
        try:
            Id_slice = search_Ids[retstart:retstart+retmax]
        except IndexError:
            Id_slice = search_Ids[retstart:]
        for Id in Id_slice: #assemble API string
            Id_fetch += f"{Id},"
        Id_fetch = Id_fetch.strip(',')
        #fetch from nuccore in INSDSeq XML format
        fetch_url=f"{efetch_base}?db=nuccore&id={Id_fetch}&rettype=gbc&rettmode=xml&api_key={api_key}"
        fetch = request_data(fetch_url)
        fetch = check_200_request_error(fetch_url, fetch, max_retries=200)
        try:
            root = ET.fromstring(fetch.text)
            entries = root.findall('./INSDSeq')
        except:
            pdb.set_trace()
        if entries is not None:
            for entry in entries:
                main, source = parse_nuccore_INSDSeq(entry)
                for key in main:
                    if key != 'acession':
                        nuccore_meta_df.loc[main['acession'], key]= main[key]
                for key in source:
                    nuccore_meta_df.loc[main['acession'], f"{key}[source]"]= source[key]
        if retstart % 60000 == 0:
            print(f"Efetch nuccore metadata {int(100*(retstart/max_count))}% complete\nSaving Progress...")
            with open(f"{output_dir}/{exp_tag}_nuccoremetadata_retstart.txt", 'w') as f:
                f.write(f'Next_nuccoremeta_retstart:{retstart+retmax}')
            nuccore_meta_df.to_csv(f"{output_dir}/{exp_tag}_nuccore_meta.csv")
    print("efetch_metatadata_nuccore complete\n")
    return nuccore_meta_df

