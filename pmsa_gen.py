""" This is the main script for generating pMSAs for Virus PPI Screening """

# Example Usage:
#   python3 scripts/pmsa_gen.py \
#       --proteome_file examples/example_data/ASFV_test_genome.fa \
#       --proteome_type raw \
#       --strain_diverse_db examples/example_data/example_db.fasta \
#       --alphafold_msa_dir exampes/example_data/af_msas/ \
#       --hhsuite_dir hhsuite/ \
#       --pmsa_output_dir examples/example_pmsa_out/
#       --cpus 20


import argparse, os, sys, pdb
import common.pmsa_pipeline as pipeline 

def main():
    parser = argparse.ArgumentParser(description="pMSA generation pipeline script.")
    parser.add_argument("--proteome_file", type=str, help="The path to \
                        the fasta formatted proteome file")
    parser.add_argument("--proteome_type", type=str, help="The type of i \
                        proteome file input. options: 'raw', 'uniprot'.")
    parser.add_argument("--strain_diverse_db", type=str, help="The path to the fasta formatted \
                         strain diverse viral protein database")
    parser.add_argument("--alphafold_msa_dir", type=str, help="The parent \
                        directory of the precomputed  AlphaFold v2.3.2 default MSA generation")
    parser.add_argument("--hhsuite_dir", type=str, help="The parent \
                        directory of the hhsuite precompiled source code")
    parser.add_argument("--pmsa_output_dir", type=str, help="The output directory \
                        of this job")
    parser.add_argument("--cpus", type=int, help="The number of CPUs available for pMSA generation")

    args = parser.parse_args()
    pipe = pipeline.pMSA_Pipeline(hhsuite_dir=args.hhsuite_dir,
                                      strain_diverse_db=args.strain_diverse_db, 
                                      proteome_file=args.proteome_file, 
                                      proteome_type=args.proteome_type,
                                      alphafold_msa_dir=args.alphafold_msa_dir, 
                                      pmsa_output_dir=args.pmsa_output_dir)            
    pipe.process(cpus=args.cpus)
    print(f"pMSA generation for {os.path.basename(args.proteome_file)} complete")

if __name__ == "__main__":
    main()
