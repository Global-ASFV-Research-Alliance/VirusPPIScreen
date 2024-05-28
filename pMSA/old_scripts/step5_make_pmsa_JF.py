import sys,os
import glob
from itertools import combinations
import pdb

def Ian_pMSA(pair):
	prot1 = pair.split('_')[0]
	prot2 = pair.split('_')[1]

	fp = open(f"{MSA_dir}/{prot1}.fas")
	sp2protsA = {}
	spsA = set([])
	isquery = 0
	for line in fp:
		if line[0] == ">":
			name = line[1:-1]
			if name == "query":
				isquery = 1
			else:
				isquery = 0
				if ":" in name:
					sp = name.split(":")[0]
				elif name[:3] == "GCF" or name[:3] == "GCA":
					sp = name.split("_")[0] + "_" + name.split("_")[1]
				elif "JGI" in name:
					sp = name.split("_")[0]
				else:
					print ("error", name)
		elif isquery:
			qseqA = line[:-1]
		else:
			seq = line[:-1]
			if seq.count("-") < len(seq) * 0.5:
				try:
					sp2protsA[sp].append([name, seq])
				except KeyError:
					sp2protsA[sp] = [[name, seq]]
					spsA.add(sp)
	fp.close()
	lenA = len(qseqA)

	fp = open(f"{MSA_dir}/{prot2}.fas")
	sp2protsB = {}
	spsB = set([])
	isquery = 0
	for line in fp:
		if line[0] == ">":
			name = line[1:-1]
			if name == "query":
				isquery = 1
			else:
				isquery = 0
				if ":" in name:
					sp = name.split(":")[0]
				elif name[:3] == "GCF" or name[:3] == "GCA":
					sp = name.split("_")[0] + "_" + name.split("_")[1]
				elif "JGI" in name:
					sp = name.split("_")[0]
				else:
					print ("error", name)
		elif isquery:
			qseqB = line[:-1]
		else:
			seq = line[:-1]
			if seq.count("-") < len(seq) * 0.5:
				try:
					sp2protsB[sp].append([name, seq])
				except KeyError:
					sp2protsB[sp] = [[name, seq]]
					spsB.add(sp)
	fp.close()
	lenB = len(qseqB)

	

	pair_count = 0
	prot1_count = 0
	prot2_count = 0
	get_sps = set([])
	rp = open(f"{pMSA_out_dir}/raw/{prot1}_{prot2}.a3m",'w')

	rp.write(">query\n")
	rp.write(qseqA + qseqB + "\n")
	for sp in spsA.intersection(spsB):
		get_sps.add(sp)
		if len(sp2protsA[sp]) == 1 and len(sp2protsB[sp]) == 1:
			pair_count += 1
			rp.write(">" + sp + "_pair\n")
			rp.write(sp2protsA[sp][0][1] + sp2protsB[sp][0][1] + "\n")


	fp.close()
	rp.close()

	os.system(f"/Users/jacobfenster/Documents/ASFV_Postdoc/hhsuite/bin/hhfilter -id 95 -cov 50 -i {pMSA_out_dir}/raw/{prot1}_{prot2}.a3m -o {pMSA_out_dir}/{prot1}_{prot2}fil.a3m")
	os.system(f"/Users/jacobfenster/Documents/ASFV_Postdoc/hhsuite/bin/hhfilter -id 80 -i {pMSA_out_dir}/{prot1}_{prot2}fil.a3m -o {pMSA_out_dir}/fil_80/{prot1}_{prot2}fil_80.a3m")
def get_fas_pairs(directory):
    # Get list of all .fas files in the directory
    fas_files = glob.glob(os.path.join(directory, "*.fas"))
    
    # Extract the basename and remove the .fas extension from each file
    basenames = [os.path.splitext(os.path.basename(f))[0] for f in fas_files]
    
    # Get all pairs of basenames
    pairs = list(combinations(basenames, 2))
    
    # Format pairs as basename1_basename2
    formatted_pairs = [f"{pair[0]}_{pair[1]}" for pair in pairs]
    
    return formatted_pairs

MSA_dir = "/Users/jacobfenster/Documents/ASFV_Postdoc/From_Ian_for_jacob/MSA_Ian"
pMSA_out_dir = "/Users/jacobfenster/Documents/ASFV_Postdoc/From_Ian_for_jacob/pMSAs"
os.makedirs(pMSA_out_dir, exist_ok=True)
os.makedirs(f"{pMSA_out_dir}/raw", exist_ok=True)
os.makedirs(f"{pMSA_out_dir}/fil_80", exist_ok=True)
pairs = get_fas_pairs(MSA_dir)
for pair in pairs:
	Ian_pMSA(pair)
