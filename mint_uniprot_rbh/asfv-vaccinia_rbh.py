# This script conducts RBH analysis between the ASFV-Georiga 2007 proteome
# and the Vaccinia WR proteome 
import os, sys, pdb
import common.rbh_tools as rbh
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar

asfv_proteome_file = "/lustrefs/fadru/projects/asfv-ppi/data/ASFV/Georgia-2007-LR743116_trans-rename.fa"
vaccinia_proteome_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
query_fasta = "/lustrefs/fadru/projects/asfv-ppi/Input/ASFV/mega_pMSA/ASFV_ASFV_hetero_homodimer/query_fastas/F334L.fasta"
db_out = "/home/jacob.fenster/scripts/rbh/tmp/" 
# going through some qc 
blastp_df = rbh._blastp_to_df(query_fasta, vaccinia_proteome_file, db_out, 
                        outfmt='10 qseqid sseqid evalue pident qcovs slen qlen qseq sseq length',
                        evalue=1)
vaccinia_proteome_db = rbh._makeblastdb(vaccinia_proteome_file, db_out, dbtype='prot')  
#blast._blastp_search(query_fasta, vaccinia_proteome_db, output_dir=db_out, 
#                       outfmt='0',
#                       sorthits=3, evalue=1, max_target_seqs=100,
#                       num_threads=1, mt_mode=1)
breakpoint()
rbh = rbh.rbh(asfv_proteome_file, vaccinia_proteome_file, db_out)
breakpoint()
