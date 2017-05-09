from scipy.stats import norm
import numpy as np

def prob_span (arg_dict, A):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	flank_len = arg_dict['flank_len']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']
	str_len = A * len(motif)
	# Compute normalization constant norm_const
	rv_dist = norm(loc = dist_mean, scale = dist_sdev)
	norm_const = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len)

	term1 = rv_dist.cdf(2 * flank_len + str_len) - rv_dist.cdf(2 * read_len + str_len)
	term2 = np.exp(-(float(2 * flank_len + str_len - dist_mean) / float(dist_sdev)) ** float(2) / float(2)) - \
			np.exp(-(float(2 * read_len + str_len - dist_mean) / float(dist_sdev)) ** float(2) / float(2))

	# CHANGE: 2*pi -> np.sqrt(2 * pi)
	comb_term = float(dist_mean - 2 * read_len - str_len) * term1 - \
				float(dist_sdev) / np.sqrt(2 * np.pi) * term2

	return comb_term / norm_const / float(2 * flank_len + str_len)





# calc_likelihood_span: compute log likelihood
def calc_likelihood_span(arg_dict, A, B, samp_r):
	dist_mean = arg_dict['read_ins_mean']
	dist_sdev = arg_dict['read_ins_stddev']
	motif = arg_dict['motif']
	ref_count = arg_dict['ref_allele_count']
	# Probability of Spanning given Allele with A copies
	prob_span__A = prob_span(arg_dict, A)
	prob_span__B = prob_span(arg_dict, B)
	prob_A = 0.5
	prob_B = 0.5

	mean_A = dist_mean - len(motif) * (A - ref_count)
	mean_B = dist_mean - len(motif) * (B - ref_count)
	rv_dist_A = norm(loc = mean_A, scale = dist_sdev)
	rv_dist_B = norm(loc = mean_B, scale = dist_sdev)

	# prob_A__span = prob_span__A * prob_A / \
	# 				(prob_span__A * prob_A + prob_span__B * prob_B)
	# prob_B__span = prob_span__B * prob_B / \
	# 				(prob_span__A * prob_A + prob_span__B * prob_B)


	# f (r | <A,B>)
	likelihood = prob_A * prob_span__A * rv_dist_A.pdf(samp_r) + \
				 prob_B * prob_span__B * rv_dist_B.pdf(samp_r)

	return likelihood
