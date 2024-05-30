""" Functions for adding viral strain diversity to AlphaFold default msa outputs """

# This script is to compute all the pMSAs for the vaccinia-vaccinia hetero homodimer set 
# need to run db searches to mirror what AF does

import os, glob, pdb, itertools
from typing import Dict, List, Tuple, Set
import pandas as pd
import numpy as np
import common.bio_parsers as bpar
import common.pmsa_tools as pmsa
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat


class pMSA_Pipeline:
    """ runs jackhmmer search of strain diverse complete viral genome database
        and combines these homologues with MSAs generated elsewhere with default AlpaFold v2.3.2 MSA pipeline 
        does all pairwise combinations of an input proteome """

    def __init__(self,
                 hhsuite_dir: str,
                 strain_diverse_db: str,
                 proteome_file: str,
                 proteome_type: str,
                 alphafold_msa_dir: str,
                 pmsa_output_dir: str
                 ):
        """ initializes the pMSA data pipeline
            Arguments:
                hhsuite_dir: path to the directory that contains all the 
                hhsuite executibles 
                strain_diverse_db: path to strain diverse database fasta 
                file 
                proteome_file: path to the viral proteome file in fasta 
                format 
                proteome_type: switch for parsing protein names from 
                proteome fasta descriptions. Optoins are "raw" and "uniprot"
                "raw" - protein name = full fasta description
                "uniprot" = protein name is P01136 for the following example 
                >sp|P01136|VGF_VACCW Pro-Viral epidermal growth factor
                alphafold_msa_dir: path to AlphaFold v2.3.2 MSA output
                parent directory
                pmsa_output_dir: path to project output directory. 
                This folder will be created at runtime

                """
        self.hhfilter_binary = f"{hhsuite_dir}bin/hhfilter"
        self.reformat_script = f"{hhsuite_dir}scripts/reformat.pl"
        self.strain_diverse_db = strain_diverse_db
        self.proteome_file = proteome_file
        self.alphafold_msa_dir = alphafold_msa_dir
        self.pmsa_output_dir = pmsa_output_dir
        self.proteome_type = proteome_type
        # do some tests to make sure the init is good

    def process(self,
                cpus: int) -> List[str]:
        """ main pMSA processing def 
            Argumnets:
                cpus: number of CPUs available to run and parallelize
                jhmmer search. Allocates two CPUs per search via threads
                and repeats the searches to fit the input cpus """ 
        
        # generate output folder structure, 
        # raises an error if the output directories already exist 
        query_fasta_dir = f"{self.pmsa_output_dir}query_fasta/"
        jhmmer_out = f"{self.pmsa_output_dir}jackhmmer_all_virus/"
        combined_msa_out = f"{self.pmsa_output_dir}combined_msas/"
        temp = f"{self.pmsa_output_dir}tmp/" 
        pmsa_out = f"{self.pmsa_output_dir}pmsas/"
        os.makedirs(self.pmsa_output_dir, exist_ok=False)
        os.makedirs(query_fasta_dir, exist_ok=False)
        os.makedirs(jhmmer_out, exist_ok=False)
        os.makedirs(combined_msa_out, exist_ok=False)
        os.makedirs(temp, exist_ok=True)
        os.makedirs(pmsa_out, exist_ok=False)
        # read proteome 
        descs, seqs = self.init_proteome()
        # generate query fasta files
        for desc, seq in zip(descs, seqs):
            with open(f"{query_fasta_dir}{desc}.fasta", "w") as f:
                f.write(f">{desc}\n{seq}")
        # search custom database with each query. Parallelize for speed
        query_fastas = glob.glob(f"{query_fasta_dir}*.fasta")
        samples = len(query_fastas)
        strain_diverse_dbs = repeat(self.strain_diverse_db, samples)
        jhmmer_outs = repeat(jhmmer_out, samples)
        temps = repeat(temp, samples)
        threads = 2
        threadss = repeat(threads, samples)
        iterations = 3 # number of jhmmer search iterations
        iterationss = repeat(iterations, samples)
        jhmmers = repeat('jackhmmer', samples)
        num_workers = cpus // threads
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            jhmmer_results = executor.map(pmsa.sub_jackhmmer_search, query_fastas, strain_diverse_dbs, jhmmer_outs, 
                                          temps, threadss, iterationss, jhmmers)
        # combine jhmmer search above with previously computed AlphaFold v2.3.2 MSA output
        protein_names = descs
        samples = len(protein_names)
        combined_msa_outs = repeat(combined_msa_out, samples)
        jhmmer_dirs = repeat(jhmmer_out, samples)
        query_fasta_dirs = repeat(query_fasta_dir, samples)
        alphafold_msa_dirs = repeat(self.alphafold_msa_dir, samples)
        temps = repeat(temp, samples)
        reformats = repeat(self.reformat_script, samples)
        num_workers = cpus // threads
        #pmsa.combined_msa(protein_names[0], combined_msa_out, jhmmer_out, query_fasta_dir, self.alphafold_msa_dir, temp, self.reformat_script)
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(pmsa.combined_msa, protein_names, combined_msa_outs, jhmmer_dirs, 
                                       query_fasta_dirs, alphafold_msa_dirs, temps, reformats)
        # generate pMSAs 

        # heterodimers
        random_pairs = list(itertools.combinations(protein_names, 2)) #all pairwise combinations of input files
        pairs = []
        # put PPIs in alphabetical order 
        for pair in random_pairs:
            new_pair = (min([pair[0], pair[1]]), max([pair[0], pair[1]]))
            pairs.append(new_pair)

        proteinA_msas, proteinB_msas = [], []
        for pair in pairs:
            proteinA_msa = f"{combined_msa_out}{pair[0]}.a3m"
            proteinB_msa = f"{combined_msa_out}{pair[1]}.a3m"
            pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99,
                              unpaired_tr=90, homodimer_tr=95, temp=temp, hhfilter=self.hhfilter_binary)
        # homodimers
        dimers = []
        proteinA_msas, proteinB_msas = [], []
        for protein in protein_names:
            dimers.append((protein, protein))
        for dimer in dimers:
            proteinA_msa = f"{combined_msa_out}{dimer[0]}.a3m"
            proteinB_msa = f"{combined_msa_out}{dimer[0]}.a3m"
            pmsa.sub_generate_pMSA(proteinA_msa, proteinB_msa, output_dir=pmsa_out, paired_tr=99,
                              unpaired_tr=90, homodimer_tr=95, temp=temp, hhfilter=self.hhfilter_binary)
        # remove the temp directory
        os.system(f"rm -r examples/example_data/example_pmsa_out")
        print(f"\npMSA generation complete. {len(pairs)} heterdimers assembled and {len(dimers)} homodimers assembled")

    def init_proteome(self):
        """ initializes the proteome descriptions and sequences """
        descs, seqs = bpar.read_fasta(self.proteome_file)
        if self.proteome_type == "raw":
            descs_out = descs
        elif self.proteome_type == "uniprot":
            desc_out = []
            for desc, seq in zip(descs, seqs):
                descs_out.append(desc.split('|')[1])
        seqs_out = []
        for seq in seqs:
            if seq[-1] == '*':
                seqs_out.append(seq[:-1])
            else:
                seqs_out.append(seqs)
        return descs_out, seqs_out

