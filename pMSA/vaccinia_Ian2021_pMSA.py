# This script takes outputs of AF and combines them with the NCBI all virus database jackhmmer search
# to generate paired and paired-unpaired pMSAs with maximal diversity 
# Process
import os, glob, pdb
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa

# Generate query fasta sequences for the vaccinia dataset 
mccraith_afmult_fastas = glob.glob(f"/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/AFmult/*.fasta")
mccraith_afmult_out_parent = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/AFmult/"
# NCBI all complete viral genomes database 
all_NCBI_db = "/lustrefs/fadru/projects/asfv-ppi/data/CompleteGenome_virus_CDS_NCBI_v2_cleanup.fasta"
# output directories  
project_dir = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia_McCraith2000_PPIs/Ian2021/"
query_fasta_dir = f"{project_dir}query_fasta/"
jhmmer_out = f"{project_dir}jackhmmer_all_virus/"
combined_msa_out = f"{project_dir}combined_msas/"
temp = f"{project_dir}tmp/"
pmsa_out = f"{project_dir}pmsas/"

os.makedirs(project_dir, exist_ok=True), os.makedirs(query_fasta_dir, exist_ok=True), os.makedirs(jhmmer_out, exist_ok=True)
os.makedirs(combined_msa_out, exist_ok=True), os.makedirs(temp, exist_ok=True)
os.makedirs(pmsa_out, exist_ok=True)
# hmmer and hh-suite software input section
jhmmer_threads = 20
reformat = "/home/jacob.fenster/scripts/pMSA/hhsuite/scripts/reformat.pl"
hhfilter = "/home/jacob.fenster/scripts/pMSA/hhsuite/bin/hhfilter"

pairs = [] # holds McCraith PPI pairs parsed from AF mult run
# extract PPI pairs and generate query fastas
if True:
    for fasta in mccraith_afmult_fastas: 
        descs, seqs = bpar.read_fasta(fasta)
        pairs.append((descs[0].split('__')[0], descs[1].split('__')[0]))
        if False:
            for seq, desc in zip(seqs, descs):
                desc = desc.split('__')[0] # get rid of duplicate proteinss that were mulitmers ie. proteinA__A/B
                with open(f"{query_fasta_dir}{desc}.fasta", 'w') as f:
                    f.write(f">{desc}\n{seq}")
# jackhmmer search of vaccinia proteins in McCraith set
# also combine msas from AF, realign, clean, and output
if False:
    for query_fasta in glob.glob(f"{query_fasta_dir}*.fasta"):
        protein_name = os.path.basename(query_fasta).split('.')[0]
        if not os.path.exists(f"{jhmmer_out}{protein_name}_jhmmer.sto"):
            jhmmer_file = pmsa.sys_jackhmmer_search(query_fasta, all_NCBI_db, jhmmer_out, threads=jhmmer_threads, iterations=3, jhmmer='jackhmmer')
        else: 
            jhmmer_file = f"{jhmmer_out}{protein_name}_jhmmer.sto"
        if not os.path.exists(f"{combined_msa_out}{protein_name}.fasta"):
            # clean up input MSAs, de-align, and concatinate 
            master_seqs, master_descs = [], [] 
            viral_descs, viral_seqs = pmsa.cleanup_sto_or_fasta(jhmmer_file, gap_ratio_tr=0.5, genome=True)
            viral_seqs = pmsa.un_align_sto_or_fasta(viral_seqs) 
            master_descs.extend(viral_descs), master_seqs.extend(viral_seqs)
            # read and combine  MSAs from AFmult output
            af_descs, af_seqs = pmsa.combine_AF_MSAs(protein_name, mccraith_afmult_out_parent, temp=temp, reformat=reformat)
            if af_descs is None:
                print(f"No AFmult msas for {protein}. skipping...")
                breakpoint()
                continue
            else:
                master_descs.extend(af_descs), master_seqs.extend(af_seqs)
                aln_descs, aln_seqs = pmsa.sys_gen_seed_hmmalign(query_fasta, master_descs, master_seqs, temp=f"{temp}aln/", phmmer='phmmer', hmmbuild='hmmbuild', hmmalign='hmmalign', reformat=reformat)
                # remove duplicate alignments prioritizing those in the all NCBI virus db
                unique_aln_descs, unique_aln_seqs = pmsa.remove_duplicates_prioitize_NCBI(aln_descs, aln_seqs)
                # send a3m MSA to output file  
                bpar.fasta_out_single_line(unique_aln_descs, unique_aln_seqs, f"{combined_msa_out}{protein_name}.fasta")
# generate McCraith pMSAs!
if True:
    for pair in pairs:
        proteinA_msa = f"{combined_msa_out}{pair[0]}.fasta"
        proteinB_msa = f"{combined_msa_out}{pair[1]}.fasta"
        # generate paired and unpaired pMSAs and ouptut to .a3m files. This returns the raw data for downstream stats
        pmsa.generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99, unpaired_tr=90, homodimer_tr=95, temp=temp,  hhfilter=hhfilter)


