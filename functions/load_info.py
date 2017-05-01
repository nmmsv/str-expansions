import pickle

def load_profile(exp_dir):
	# loadinf profile
	try:
		with open(exp_dir + '/profile.txt', 'r') as f:
			arg_dict = pickle.load(f)
	except:
		raise
	return arg_dict

def extract_locus_info(locus_bed):
	with open(locus_bed, 'r') as f:
		row = f.readline().split()
		chrom = row[0]
		locus_start = int(row[1])
		locus_end = int(row[2])
	return chrom, locus_start, locus_end