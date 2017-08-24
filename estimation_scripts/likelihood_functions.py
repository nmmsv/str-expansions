from scipy.stats import norm
import numpy as np
import sys


sys.path.append('/storage/nmmsv/str-expansions/functions/')
from load_info import load_profile, extract_locus_info
##############################################
############# FRR READ CLASS #################
def FRR_class_prob_old (arg_dict, A):
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
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

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

def FRR_class_prob (arg_dict, A):
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
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

	coef0 = 1.0 / norm_const / (2.0 * flank_len + str_len - 2.0 * read_len)
	coef1 = - float(dist_sdev ** 2)
	term1 = rv_dist.pdf(str_len) - rv_dist.pdf(2 * read_len)
	coef2 = dist_mean - read_len
	term2 = rv_dist.cdf(str_len) - rv_dist.cdf(2 * read_len)
	coef3 = str_len - read_len
	term3 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(str_len)
	
	if str_len >= 2 * read_len:
		return_value = coef0 * (coef1 * term1 + coef2 * term2 + coef3 * term3)
	else:
		return_value = coef0 * coef3 * term3
	# print A, return_value
	return return_value

def FRR_read_prob (arg_dict, A, sample_dfl):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	flank_len = arg_dict['flank_len']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']

	str_len = A * len(motif)
	# Commenting lines below cause redundant
	# if str_len < read_len:		# condition: L > r for this read to be possible
	# 	return 0
	# insert size distribution
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)

	# Compute normalization constant norm_const
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

	term1 = rv_dist.cdf(read_len + sample_dfl + str_len) - rv_dist.cdf(2 * read_len + sample_dfl)

	return_value = 1 / norm_const * term1

	return return_value


def FRR_allele_likelihood (arg_dict, A, sample_ins):
	return FRR_class_prob(arg_dict, A) * FRR_read_prob(arg_dict, A, sample_ins)
def FRR_genotype_likelihood(arg_dict, A, B, sample_ins): # this is not used -> implementation embedded in sam ll functions
	prob_A = 0.5
	prob_B = 0.5
	return prob_A * FRR_allele_likelihood(arg_dict, A, sample_ins) + \
			prob_B * FRR_allele_likelihood(arg_dict, B, sample_ins)
##############################################
########### SPANNING READ CLASS ##############
def span_class_prob (arg_dict, A):
	dist_mean = arg_dict['read_ins_mean'] 		# (\mu)
	dist_sdev = arg_dict['read_ins_stddev']		# (\sigma)
	flank_len = arg_dict['flank_len']			# (F)
	read_len = arg_dict['read_len']				# (r)
	motif = arg_dict['motif']
	str_len = A * len(motif)					# (L)
	# Compute normalization constant norm_const (C)
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

	coef0 = 1.0 / norm_const / float(2.0 * flank_len + str_len - 2.0 * read_len)

	coef1 = float(dist_mean - str_len)
	term1 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(max(2 * read_len, str_len))
	coef2 = - float(dist_sdev ** 2)
	term2 = rv_dist.pdf(2 * flank_len + str_len) - rv_dist.pdf(max(2 * read_len, str_len))
	# print 
	# print coef2
	return_value = coef0 * (coef1 * term1 + coef2 * term2)

	# print return_value, coef0, coef1, coef2
	# print return_value, term1, term2
	# print
	## OLD VERSION!
	# coef0 = 1.0 / norm_const / float(2 * flank_len + str_len)
	# coef1 = - float(dist_sdev)  / np.sqrt(2 * np.pi)
	# term1 = np.exp(-0.5 * (float(2 * flank_len + str_len - dist_mean) / float(dist_sdev)) ** float(2)) - \
	# 		np.exp(-0.5 * (float(2 * read_len + str_len - dist_mean) / float(dist_sdev)) ** float(2))
	# coef2 = float(dist_mean - 2 * read_len - str_len)
	# term2 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len + str_len)
	# return_value = coef0 * (coef1 * term1 + coef2 * term2)
	
	# print coef1, term1, coef1 * term1, '\t', coef2, term2, coef2 * term2
	return return_value


def span_read_prob (arg_dict, A, sample_ins):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	motif = arg_dict['motif']
	ref_count = arg_dict['ref_allele_count']

	mean_A = dist_mean - len(motif) * (A - ref_count)
	# shifted insert size distribution
	# NOTE NOTE NOTE DEVIDIED BY 3!!
	rv_dist_A = norm(loc = mean_A, scale = dist_sdev)
	# print np.log(rv_dist_A.pdf(sample_ins))
	return rv_dist_A.pdf(sample_ins)

def span_allele_likelihood (arg_dict, A, sample_ins):  
	# print np.log(span_class_prob(arg_dict, A) * span_read_prob(arg_dict, A, sample_ins))
	return span_class_prob(arg_dict, A) * span_read_prob(arg_dict, A, sample_ins)
def span_genotype_likelihood(arg_dict, A, B, sample_ins): # this is not used -> implementation embedded in sam ll functions
	prob_A = 0.5
	prob_B = 0.5
	return prob_A * span_allele_likelihood(arg_dict, A, sample_ins) + \
			prob_B * span_allele_likelihood(arg_dict, B, sample_ins)
##############################################
##############################################
########### ENCLOSING READ CLASS #############
def encl_class_prob(arg_dict, A):
	flank_len = arg_dict['flank_len']			# (F)
	read_len = arg_dict['read_len']				# (r)
	motif = arg_dict['motif']
	str_len = A * len(motif)					# (L)

	if read_len < str_len:
		return 0.0

	return_value = float(read_len - str_len) / float(2 * flank_len + str_len - 2 * read_len)
	return return_value

def encl_read_prob(arg_dict,A, sample_ncopy):
	# FILLLLLL INNNNN Stutter model
	u = 0.01
	d = 0.02
	p = 0.95
	delta = sample_ncopy - A
	if delta == 0:
		return 1 - u - d
	elif delta > 0:
		return u * p * (1.0 - p) ** (delta - 1.0)
	else:
		return d * p * (1.0 - p) ** (-delta - 1.0)

	# # For now, I'll use normal distibution:
	# encl_sdev = 1.0
	# rv_encl = norm(loc = sample_ncopy, scale = encl_sdev)
	# return rv_encl.pdf(A)

def encl_allele_likelihood(arg_dict, A, sample_ncopy):
	return encl_class_prob(arg_dict, A) * encl_read_prob(arg_dict, A, sample_ncopy)

def encl_genotype_likelihood(arg_dict, A, B, sample_ins): # this is not used -> implementation embedded in sam ll functions
	prob_A = 0.5
	prob_B = 0.5
	return prob_A * encl_allele_likelihood(arg_dict, A, sample_ins) + \
			prob_B * encl_allele_likelihood(arg_dict, B, sample_ins)

def test_values_for_C(arg_dict, allele, sample_ins, sample_dfl, sample_ncopy):
	print 'allele:\t\t\t', allele
	print 'dist_mean:\t\t', arg_dict['read_ins_mean']
	print 'dist_sdev:\t\t', arg_dict['read_ins_stddev']
	print 'motif_len:\t\t', len(arg_dict['motif'])
	print 'flank_len:\t\t', arg_dict['flank_len']
	print 'read_len: \t\t', arg_dict['read_len']
	print 'ref_count:\t\t', arg_dict['ref_allele_count']
	print '########## SPAN ##########'
	print 'data(sample_ins):\t', sample_ins
	print 'Class prob:\t\t', np.log(span_class_prob(arg_dict, allele))
	print 'Read prob:\t\t', np.log(span_read_prob(arg_dict, allele, sample_ins))
	print 'Allele ll:\t\t', np.log(span_allele_likelihood(arg_dict, allele, sample_ins))
	print '########## RFF ###########'
	print 'data(sample_dfl):\t', sample_dfl
	print 'Class prob:\t\t', np.log(FRR_class_prob(arg_dict, allele))
	print 'Read prob:\t\t', np.log(FRR_read_prob(arg_dict, allele, sample_dfl))
	print 'Allele ll:\t\t', np.log(FRR_allele_likelihood(arg_dict, allele, sample_dfl))
	print '########## ENCL ##########'
	print 'data(sample_ncopy):\t', sample_ncopy
	print 'Class prob:\t\t', np.log(encl_class_prob(arg_dict, allele))
	print 'Read prob:\t\t', np.log(encl_read_prob(arg_dict, allele, sample_ncopy))
	print 'Allele ll:\t\t', np.log(encl_allele_likelihood(arg_dict, allele, sample_ncopy))

def test_sum_ll_for_C(arg_dict, allele1, allele2, data_v):
	ll_sum_encl = 0;
	ll_sum_span = 0;
	ll_sum_frr = 0;
	for datum in data_v:
		ll_sum_encl += np.log(encl_genotype_likelihood(arg_dict, allele1, allele2, datum))
		ll_sum_span += np.log(span_genotype_likelihood(arg_dict, allele1, allele2, datum))
		ll_sum_frr += np.log(FRR_genotype_likelihood(arg_dict, allele1, allele2, datum))
	print 'Enclosing sum_ll:\t\t', ll_sum_encl
	print 'Spanning sum_ll:\t\t', ll_sum_span
	print 'FRR sum_ll:\t\t\t', ll_sum_frr
##############################################

exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_18_class2_cov50_dist400/'
arg_dict = load_profile(exp_dir)

# A = 30
# B = 50
# sample_ins = 300

# print span_genotype_likelihood(arg_dict, A, B, sample_ins)
# print FRR_genotype_likelihood(arg_dict, A, B, 400)

# for i in range(10,200,50):
# 	for j in range(10,200,50):
# 		print i, j, FRR_genotype_likelihood(arg_dict, i, j, 100)

# test_values_for_C(arg_dict, 55, 430, 80, 24)
# test_sum_ll_for_C(arg_dict, 20, 50, [10,20,30,40])
weights = {	'frr': 1.0, \
			'srp': 1.0, \
			'er': 1.0, \
			'allele_1': 0.5,\
			'allele_2': 0.5}
x = [10, 30]
print (-weights['er'] * np.log(encl_genotype_likelihood(arg_dict, x[0], x[1], 10)) + \
					-weights['srp'] * np.log(span_genotype_likelihood(arg_dict, x[0], x[1], 380)))
x = [45, 65]
print (-weights['srp'] * np.log(span_genotype_likelihood(arg_dict, x[0], x[1], 380)) + \
					-weights['frr' ] * np.log(FRR_genotype_likelihood(arg_dict, x[0], x[1], 80)))
