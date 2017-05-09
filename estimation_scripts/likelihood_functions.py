from scipy.stats import norm
import numpy as np
import sys


sys.path.append('/storage/nmmsv/str-expansions/functions/')
from load_info import load_profile, extract_locus_info
##############################################
############# IRR READ CLASS #################
def IRR_class_prob (arg_dict, A):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	flank_len = arg_dict['flank_len']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']
	str_len = A * len(motif)

	if str_len < read_len:		# condition: L > r for this read to be possible
		return 0
	# insert size distribution
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)

	# Compute normalization constant norm_const
	norm_const = rv_dist.cdf(2 * flank_len - str_len) - rv_dist.cdf(2 * read_len)

	coef0 = 1 / norm_const / (2 * flank_len + str_len)
	coef1 = - dist_sdev / np.sqrt(2 * np.pi)
	term1 = np.exp(-((str_len + read_len - dist_mean) / dist_sdev) ** 2 / 2) - \
			np.exp(-((2 * read_len - dist_mean) / dist_sdev) ** 2 / 2)
	coef2 = dist_mean - 2 * read_len
	term2 = rv_dist.cdf(str_len + read_len) - rv_dist.cdf(2 * read_len)
	coef3 = str_len - read_len
	term3 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(str_len + read_len)
	
	return_value = coef0 * (coef1 * term1 + coef2 * term2 + coef3 * term3)

	return return_value


def IRR_read_prob (arg_dict, A, sample_dfl):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	flank_len = arg_dict['flank_len']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']

	str_len = A * len(motif)

	if str_len < read_len:		# condition: L > r for this read to be possible
		return 0
	# insert size distribution
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)

	# Compute normalization constant norm_const
	norm_const = rv_dist.cdf(2 * flank_len - str_len) - rv_dist.cdf(2 * read_len)

	term1 = rv_dist.cdf(read_len + sample_dfl + str_len) - rv_dist.cdf(2 * read_len + sample_dfl)

	return_value = 1 / norm_const * term1

	return return_value


def IRR_allele_likelihood (arg_dict, A, sample_ins):
	return IRR_class_prob(arg_dict, A) * IRR_read_prob(arg_dict, A, sample_ins)
def IRR_genotype_likelihood(arg_dict, A, B, sample_ins):
	prob_A = 0.5
	prob_B = 0.5
	return prob_A * IRR_allele_likelihood(arg_dict, A, sample_ins) + \
			prob_B * IRR_allele_likelihood(arg_dict, B, sample_ins)
##############################################
########### SPANNING READ CLASS ##############
def span_class_prob (arg_dict, A):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	flank_len = arg_dict['flank_len']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']
	str_len = A * len(motif)
	# Compute normalization constant norm_const
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

	coef0 = 1 / norm_const / float(2 * flank_len + str_len)

	coef1 = - float(dist_sdev)  / np.sqrt(2 * np.pi)
	term1 = np.exp(-0.5 * (float(2 * flank_len + str_len - dist_mean) / float(dist_sdev)) ** float(2)) - \
			np.exp(-0.5 * (float(2 * read_len + str_len - dist_mean) / float(dist_sdev)) ** float(2))
	coef2 = float(dist_mean - 2 * read_len - str_len)
	term2 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len + str_len)

	return_value = coef0 * (coef1 * term1 + coef2 * term2)

	return return_value


def span_read_prob (arg_dict, A, sample_ins):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	motif = arg_dict['motif']
	ref_count = arg_dict['ref_allele_count']

	mean_A = dist_mean - len(motif) * (A - ref_count)
	# shifted insert size distribution
	rv_dist_A = norm(loc = mean_A, scale = dist_sdev)

	return rv_dist_A.pdf(sample_ins)

def span_allele_likelihood (arg_dict, A, sample_ins):
	return span_class_prob(arg_dict, A) * span_read_prob(arg_dict, A, sample_ins)
def span_genotype_likelihood(arg_dict, A, B, sample_ins):
	prob_A = 0.5
	prob_B = 0.5
	return prob_A * span_allele_likelihood(arg_dict, A, sample_ins) + \
			prob_B * span_allele_likelihood(arg_dict, B, sample_ins)
##############################################
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_18_class2_cov50_dist400/'
arg_dict = load_profile(exp_dir)

# A = 30
# B = 50
# sample_ins = 300

# print span_genotype_likelihood(arg_dict, A, B, sample_ins)
# print IRR_genotype_likelihood(arg_dict, A, B, 400)

# for i in range(10,200,50):
# 	for j in range(10,200,50):
# 		print i, j, IRR_genotype_likelihood(arg_dict, i, j, 100)
