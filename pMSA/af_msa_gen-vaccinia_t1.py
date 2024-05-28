# This is to run the msa generation only alphafold data pipeline via singularity
import os, glob, pdb, time, itertools, subprocess, sys
import concurrent, multiprocessing
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar
import pmsa_tools as pmsa



def generate_af_msas(fasta_path, 
                     output_dir, 
                     container='/lustrefs/software/public/containers/alphafold_2.3.2-1.sif', 
                     msa_script='/lustrefs/software/public/apps/alphafold-2.3.2/alphafold-2.3.2/run_alphafold_msa_only.py'):
    print(f"msa run for {fasta_path} started")
    timei = time.time()
    binds = f"{msa_script}:/app/alphafold/run_alphafold_msa_only.py"
    commands = f"python3 run_alphafold_msa_only.py \
--fasta_paths={fasta_path} \
--output_dir={output_dir} \
--data_dir=/lustrefs/public/databases/alphafold/ \
--model_preset=multimer \
--uniref90_database_path=/lustrefs/public/databases/alphafold/uniref90/uniref90.fasta \
--mgnify_database_path=/lustrefs/public/databases/alphafold/mgnify/mgy_clusters_2022_05.fa \
--template_mmcif_dir=/lustrefs/public/databases/alphafold/pdb_mmcif/mmcif_files/ \
--max_template_date=2020-01-01 \
--obsolete_pdbs_path=/lustrefs/public/databases/alphafold/pdb_mmcif/obsolete.dat \
--use_gpu_relax=False \
--db_preset=full_dbs \
--bfd_database_path=/lustrefs/public/databases/alphafold/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
--uniref30_database_path=/lustrefs/public/databases/alphafold/uniref30/UniRef30_2021_03 \
--pdb_seqres_database_path=/lustrefs/public/databases/alphafold/pdb_seqres/pdb_seqres.txt \
--uniprot_database_path=/lustrefs/public/databases/alphafold/uniprot/uniprot.fasta"
    singularity_command = f"singularity exec --pwd /app/alphafold --bind {binds} {container} {commands}"
    msa_gen = subprocess.run(singularity_command, shell=True, text=True, capture_output=True)
    if msa_gen.returncode != 0:
        print(f"ERROR with {fasta_path}, returncode={msa_gen.returncode} : {msa_gen.stderr}")
        print('')
        return
    timef = time.time()
    print(f"msa run complete for {fasta_path}")
    print(f"Runtime: {(timef-timei)/60:.1f} minutes")
    print('')
    


if __name__ == '__main__':
    # generate fasta homodimer alphafold multimer input fastas from uniprot vaccinia WR proteome 
    vaccinia_proteome_file = "/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/Vaccinia_Virus-WR_uniprot_ref_proteome.fasta"
    fasta_input = "/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/homodimer_msa_gen/tranche1/"
    descs, seqs = bpar.read_fasta(vaccinia_proteome_file)
    if False:
        for desc, seq in zip(descs, seqs):
            protein = f"{desc.split('|')[1]}.{desc.split('|')[2].split()[0]}"
            basename = desc.split('|')[1]
            with open(f"{fasta_input}{basename}__{basename}.fasta", 'w') as f:
                f.write(f">{protein}__A\n{seq}\n>{protein}__B\n{seq}")
    if True:
        output_dir = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/homodimer/AFmult_msa"
        fasta_paths = glob.glob(f"{fasta_input}*.fasta")
        output_dirs = repeat(output_dir, len(fasta_paths))
        num_workers = 12
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            results = executor.map(generate_af_msas, fasta_paths, output_dirs)
