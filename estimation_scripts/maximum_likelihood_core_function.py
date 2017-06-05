from likelihood_spanning import calc_likelihood_span
import csv, sys
import numpy as np
from scipy.optimize import minimize
import pickle, argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile, extract_locus_info
from likelihood_functions import span_genotype_likelihood, FRR_genotype_likelihood
from likelihood_functions import span_allele_likelihood, FRR_allele_likelihood
from realignment import expansion_aware_realign
from extract_genome import extract_pre_post_flank

def extract_col_sam(columns, tag):
	for item in columns:
		parts = item.split(':')
		if len(parts) != 3:
			pass
		elif parts[0] == tag:
			if parts[1] == 'i':
				return int(parts[2])
			else:
				return parts[2]

	raise ValueError ('Cannot find tag: ' + tag)

def FRR_sam_likelihood (A, sam_path, arg_dict):
	locus = arg_dict['locus']
	chrom, locus_start, locus_end = extract_locus_info(locus)
	log_likelihood = 1
	with open(sam_path, 'r') as irr_handle:
		for record in csv.reader(irr_handle, dialect = 'excel-tab'):
			if record[0][0] != '@':
				sample_dfl = extract_col_sam(record, 'om')
				samp_likelihood = FRR_allele_likelihood(arg_dict, A, sample_dfl)
				# print record[0], '\t', locus_start - int(record[3]), '\t', samp_likelihood
				if samp_likelihood > 0:
					samp_log_likelihood = np.log(samp_likelihood)
				else:
					samp_log_likelihood = -50
				log_likelihood = log_likelihood + samp_log_likelihood
				# log_likelihood = log_likelihood * samp_likelihood
	return log_likelihood


def span_sam_likelihood (A, sam_path, arg_dict):
	locus = arg_dict['locus']
	mean_ins_size = arg_dict['read_ins_mean']
	chrom, locus_start, locus_end = extract_locus_info(locus)
	log_likelihood = 0
	yo = 0
	nn = 0
	with open(sam_path, 'r') as irr_handle:
		for record in csv.reader(irr_handle, dialect = 'excel-tab'):
			if record[0][0] != '@' and int(record[8]) != 0:
				sample_ins = extract_col_sam(record, 'is')
				samp_likelihood = span_allele_likelihood(arg_dict, A, sample_ins)

				# print record[0], '\t', sample_ins, '\t', np.abs(sample_ins - mean_ins_size) / 3 + 10, '\t', samp_likelihood
				yo = yo + np.abs(sample_ins - mean_ins_size) / 3 + 10
				nn = nn + 1
				if samp_likelihood > 0:
					samp_log_likelihood = np.log(samp_likelihood)
				elif np.abs(samp_likelihood) < 10**-20:		# accounting for comutational errors
					samp_log_likelihood = np.log(np.abs(samp_likelihood))
				elif samp_likelihood == 0:
					samp_log_likelihood = -50
				else:
					print 'Error! Negative likelihood:', samp_likelihood
				log_likelihood = log_likelihood + samp_log_likelihood
	# print yo / nn
	return log_likelihood

def encl_sam_genotype (sam_path, arg_dict):
	exp_dir = arg_dict['exp_dir']
	read_len = arg_dict['read_len']

	nCopy_dict = {}
	total_count = 0
	with open(sam_path, 'r') as encl_handle:
		for record in csv.reader(encl_handle, dialect = 'excel-tab'):
			if record[0][0] != '@':
				sample = record[9]
				nCopy = extract_col_sam(record, 'nc')
				if nCopy not in nCopy_dict:
					nCopy_dict[nCopy] = 1
				else:
					nCopy_dict[nCopy] = nCopy_dict[nCopy] + 1
				total_count = total_count + 1
			nCopy_list = nCopy_dict.keys()
			freq_list = []
			for nCopy in nCopy_list:
				freq_list.append(float(nCopy_dict[nCopy]) / float(total_count))

	return nCopy_list, freq_list


def ml_enclosing_spanning(in_pref, exp_dir):
	arg_dict = load_profile(exp_dir)

	in_span = in_pref + '_srp.sam'
	in_encl = in_pref + '_er.sam'
	in_frep = in_pref + '_frr.sam'


	fn = lambda x: (-1 * FRR_sam_likelihood(x[0], in_frep, arg_dict) + \
					-1 * span_sam_likelihood(x[0], in_span, arg_dict))

	fn_span = lambda x: (-1 * span_sam_likelihood(x[0], in_span, arg_dict))

	fn_frep = lambda x: (-1 * FRR_sam_likelihood(x[0], in_frep, arg_dict))

	
	# def diminisher(x, encl_allele):
	# 	val = fn_span_single(encl_allele)
	# 	sig = 1
	# 	amplitude_damper = 0.01
	# 	return amplitude_damper * val * np.exp(-0.5 * ((x - encl_allele[0]) / sig) ** 2)
	# fn_span_single_diminished = lambda x: fn_span_single(x) + diminisher(x[0], encl_allele)
	# Now try using spanning read pairs:
	res = minimize(fn_span, x0 = [30], \
					method = 'L-BFGS-B', \
					bounds = [(1, 250)])
	return [int(round(j)) for j in res.x] 

def ml_FRR_spanning(in_pref, exp_dir, weights):
	arg_dict = load_profile(exp_dir)
	in_span = in_pref + '_srp.sam'
	in_encl = in_pref + '_er.sam'
	in_frep = in_pref + '_frr.sam'
	nc = 80
	fn = lambda x: (-weights['frr'] * FRR_sam_likelihood(x[0], in_frep, arg_dict) + \
					-weights['srp'] * span_sam_likelihood(x[0], in_span, arg_dict))
	res = minimize(fn, x0 = [70], \
						method = 'L-BFGS-B', \
						bounds = [(1, 250)])
	return [int(round(j)) for j in res.x] 



# print ml_enclosing_spanning('nc_120', '../')

# arg_dict = load_profile('../')
# nc = 80
# fn_frep = lambda x: (-1 * FRR_sam_likelihood(x[0], 'nc_'+str(nc)+'_frr.sam', arg_dict))
# fn = lambda x: (-0.1 * FRR_sam_likelihood(x[0], 'nc_'+str(nc)+'_frr.sam', arg_dict) + \
# 				-1 * span_sam_likelihood(x[0], 'nc_'+str(nc)+'_srp.sam', arg_dict))
# res = minimize(fn, x0 = [70], \
# 					method = 'L-BFGS-B', \
# 					bounds = [(1, 250)])
# print [int(round(j)) for j in res.x] 
# for i in range (5,50, 5):
# 	print i, '\t', span_sam_likelihood(i, 'nc_40_srp.sam', arg_dict)

# print span_sam_likelihood(70, 'nc_70_srp.sam', arg_dict)

