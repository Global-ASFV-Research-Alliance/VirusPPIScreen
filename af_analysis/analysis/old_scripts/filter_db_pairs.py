import os, sys, shutil
import pandas as pd

def filter_files_by_protein_pairs(file_list, protein_pairs):
    filtered_files = []
    not_found_protein_pairs = set(protein_pairs)
    for file in file_list:
        for proteinA, proteinB in protein_pairs:
            if proteinA in file and proteinB in file:
                filtered_files.append(file)
                if (proteinA, proteinB) in not_found_protein_pairs:
                    not_found_protein_pairs.remove((proteinA, proteinB))
    return filtered_files, not_found_protein_pairs

def copy_files(source_dir, dest_dir, files_to_copy):
    for file in files_to_copy:
        shutil.copy(os.path.join(source_dir, file), os.path.join(dest_dir, file))

def copy_pair_files_by_filename(pMSA_dir, protein_pairs_db, dest_dir):
    """
    This def pulls pMSA files and copies them to dest_dir given the pMSA dir, the path to the paired dir,
    copies pMSAs to dest dir
    """
    pMSA_filelist = os.listdir(pMSA_dir)
    os.makedirs(dest_dir, exist_ok=True)
    pair_db = pd.read_csv(protein_pairs_db)
    protein_pairs = [(row['proteinA'], row['proteinB']) for _, row in pair_db.iterrows()]
    filtered_pMSAs, not_found = filter_files_by_protein_pairs(pMSA_filelist, protein_pairs)
    print(f"There are {int(100*float(len(filtered_pMSAs))/float(len(protein_pairs)))}% files returned from the given directory")
    if not_found:
        print("Protein pairs not found in any files:")
        for geneA, geneB in not_found:
            print(f"{geneA}, {geneB}")
    copy_files(pMSA_dir, dest_dir, filtered_pMSAs)

def filter_data_by_pairs(data_file, protein_pairs_db, dest_dir):
    data = pd.read_csv(data_file, index_col=0)
    filtered_data = pd.DataFrame(columns=data.columns)
    pair_db = pd.read_csv(protein_pairs_db)
    protein_pairs = [(row['proteinA'], row['proteinB']) for _, row in pair_db.iterrows()]
    not_found_protein_pairs = set(protein_pairs)
    for index in data.index:
        for proteinA, proteinB in protein_pairs:
            if proteinA in index and proteinB in index:
                filtered_data.loc[index] = data.loc[index]
                if (proteinA, proteinB) in not_found_protein_pairs:
                    not_found_protein_pairs.remove((proteinA, proteinB))
    print(f"There are {int(100*float(len(filtered_data.index))/float(len(protein_pairs)))}% files returned from the given directory")
    if not_found_protein_pairs:
        print("Protein pairs not found in any index of data:")
        for geneA, geneB in not_found_protein_pairs:
            print(f"{geneA}, {geneB}")
    filtered_data.to_csv(f"{dest_dir}/{os.path.basename(data_file).split('.')[0]}_fil_{os.path.basename(protein_pairs_db).split('.')[0]}.csv")
    return filtered_data, not_found_protein_pairs
#usage 1: pytyon3 filter_db_pairs.py copy_pairs /path/to/pMSA/files /path/to/pairdb.csv /path/to/output/
#usage 2: python3 filter_db_pairs.py filter_db /path/to/data.csv /path/to/pairdb.csv /path/to/output/
if sys.argv[1] == 'copy_pairs':
    pMSA_dir, protein_pairs_db, dest_dir = sys.argv[2], sys.argv[3], sys.argv[4]
    copy_pair_files_by_filename(pMSA_dir, protein_pairs_db, dest_dir)
elif sys.argv[1] == 'filter_db':
    data_file, protein_pairs_db, dest_dir = sys.argv[2], sys.argv[3], sys.argv[4]
    filtered_data, not_found_protein_pairs = filter_data_by_pairs(data_file, protein_pairs_db, dest_dir)
else:
    print('invalid switch input.\nusage: python3 analyze_db_pairs filter_db /path/to/data.csv /path/to/pairdb.csv /path/to/output/')

