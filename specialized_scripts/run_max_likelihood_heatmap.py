
import subprocess
import sys

sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from maximum_likelihood_core_function import ml_enclosing_spanning
from load_info import load_profile, extract_locus_info

import argparse
parser = argparse.ArgumentParser('Estimate STR length, calculate error, store values for heatmap draw script')
parser.add_argument('--exp-dir',	type = str, required = True)
# parser.add_argument('--out-path', 	type = str, required = True)


args = parser.parse_args()

exp_dir = args.exp_dir
# out_path = args.out_path

algn_read_dir = exp_dir + '/aligned_read/'

arg_dict = load_profile(exp_dir)

copy_list = arg_dict['num_copy']
cov_list = arg_dict['coverage']
repo_dir = arg_dict['repo_dir']

ml_script_dir = repo_dir + '/estimation_scripts/Maximum_Likelihood.py'


for coverage in cov_list:
	for nc in copy_list:
		in_pref = algn_read_dir + 'nc' + str(nc) + '_cv' + str(coverage)
		genot = ml_enclosing_spanning(in_pref, exp_dir)
		print 'cov:', coverage, '\t nCopy:', nc, '\tGenotype:', ml_enclosing_spanning(in_pref, exp_dir)
		alt_allele = genot[0]
		

print copy_list
print cov_list