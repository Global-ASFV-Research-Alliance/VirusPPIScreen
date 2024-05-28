import os, sys, glob, pdb
import pandas as pd
import analysis.AFppi as ppi
sys.path.append("/home/jacob.fenster/scripts/pMSA/common/")
import bio_parsers as bpar
# input data
mccraith_afmult = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/AFmult/"
mccraith_ian_paired = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/Ian2021/paired/"
mccraith_ian_paired_unpaired = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/Ian2021/paired_unpaired/"
mccraith_ian_homodimer = "/lustrefs/fadru/projects/asfv-ppi/Output/Vaccinia/McCraith2000/Ian2021/homodimer/"

vacc_out = "/lustrefs/fadru/projects/asfv-ppi/Output/AF_Virus_analysis/Vaccinia/"
# Generate McCraith 2000 protein length df from fasta inputs
mccraith_afmult_fastas = glob.glob("/lustrefs/fadru/projects/asfv-ppi/Input/Vaccinia/AFmult/*.fasta")
mccraith_protlen_df = pd.DataFrame(columns=['length'])
for fasta in mccraith_afmult_fastas:
    descs, seqs = bpar.read_fasta(fasta)
    protA, protB, protAlen, protBlen = descs[0].split('__')[0], descs[1].split('__')[0], len(seqs[0]), len(seqs[1])
    if protA not in mccraith_protlen_df.index:
        mccraith_protlen_df.loc[protA, 'length'] = protAlen
    if protB not in mccraith_protlen_df.index:
        mccraith_protlen_df.loc[protB, 'length'] = protBlen
mccraith_protlen_df.to_csv("/lustrefs/fadru/projects/asfv-ppi/data/Vaccinia/McCraith2000_proteinlen.csv")

mccraith_afmult_df = ppi.AFmult_data_extract(mccraith_afmult)
mc_i_p_df  = ppi.dir_max_contact(mccraith_ian_paired, mccraith_protlen_df)
mc_i_pu_df  = ppi.dir_max_contact(mccraith_ian_paired_unpaired, mccraith_protlen_df)
mc_i_ho_df = ppi.dir_max_contact(mccraith_ian_homodimer, mccraith_protlen_df)
mccraith_afmult_df.to_csv(f"{vacc_out}McCraith2000_AFmult.csv")
mc_i_p_df.to_csv(f"{vacc_out}McCraith2000_Ian_heterodimer_paired.csv")
mc_i_pu_df.to_csv(f"{vacc_out}McCraith2000_Ian_heterodimer_paired-unpaired.csv")
mc_i_ho_df.to_csv(f"{vacc_out}McCraith2000_Ian_homodimer.csv")

