# Steps:
# First run genotype enclosing -> One or two estimated alleles
# Then run ML in hipstr area -> Two estimated alleles
# Use the result from previous steps for 1-D maximum likelihood -> Up to three estimated alleles
# Perform ML in non-hipstr region -> Two estimated alleles
# (We have up to 9 alleles) -> delete duplicates (keep doubles)
# Check all combinations to find maximum likelihood

import math
import random, sys
from simanneal import Annealer
import numpy as np
from scipy.optimize import minimize

sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile
from maximum_likelihood_core_function import span_sam_likelihood, FRR_sam_likelihood, encl_sam_likelihood, ml_encl_FRR_span, encl_sam_genotype

def top_two(allele_list, freq_list):
	max_freq_1 = 0
	max_freq_2 = 0
	max_allele_1 = 0
	max_allele_2 = 0
	j = 0
	for freq in freq_list:
		if freq > max_freq_1:
			max_freq_2 = max_freq_1
			max_freq_1 = freq
			max_allele_2 = max_allele_1
			max_allele_1 = allele_list[j]
		elif freq > max_freq_2:
			max_freq_2 = freq
			max_allele_2 = allele_list[j]
		j = j + 1

	return_value = [max_allele_1, max_allele_2]
	if return_value == [0, 0]:
		print 'Weird [0,0] encountered in top_two() function.'

	return return_value

class OptimalAlleleLengthProblem(Annealer):
	def __init__(self, state, in_pref, weights, exp_dir):
		self.random_range = 5
		self.in_span = in_pref + '_srp.sam'
		self.in_encl = in_pref + '_er.sam'
		self.in_frep = in_pref + '_frr.sam'
		self.weights = weights
		self.arg_dict = load_profile(exp_dir)
		super(OptimalAlleleLengthProblem, self).__init__(state)  # important!

	def move(self):
		change0 = random.randint(-self.random_range, self.random_range)
		change1 = random.randint(-self.random_range, self.random_range)
		# print self.state
		search_range = self.arg_dict['read_len'] / len(self.arg_dict['motif']) + 2
		self.state[0] = (self.state[0] + change0) % search_range
		self.state[1] = (self.state[1] + change1) % search_range

	def energy(self):
		e = (-self.weights['frr'] * FRR_sam_likelihood (self.state, self.in_frep, self.arg_dict, self.weights) + \
						-self.weights['srp'] * span_sam_likelihood(self.state, self.in_span, self.arg_dict, self.weights) + \
						-self.weights['er' ] * encl_sam_likelihood(self.state, self.in_encl, self.arg_dict, self.weights))
		return e

def run_sim_annealing(exp_dir, in_pref, weights, oalp_params):
	arg_dict = load_profile(exp_dir)
	oalp = OptimalAlleleLengthProblem(oalp_params['init_state'], in_pref, weights, exp_dir)
	oalp.steps = oalp_params['steps']
	oalp.updates = oalp_params['updates']
	oalp.Tmax = oalp_params['Tmax']
	oalp.Tmin = oalp_params['Tmin']
	oalp.copy_strategy = oalp_params['copy_strategy']
	state, e = oalp.anneal()
	return state, e

def run_max_likelihood(in_pref, exp_dir, weights, ml_params):
	arg_dict = load_profile(exp_dir)
	in_span = in_pref + '_srp.sam'
	in_encl = in_pref + '_er.sam'
	in_frep = in_pref + '_frr.sam'
	if ml_params['dimensions'] == 1:
		fixed_allele = ml_params['fixed_allele']
		fn = lambda x: (-weights['frr'] * FRR_sam_likelihood ([fixed_allele, x[0]], in_frep, arg_dict, weights) + \
						-weights['srp'] * span_sam_likelihood([fixed_allele, x[0]], in_span, arg_dict, weights) + \
						-weights['er' ] * encl_sam_likelihood([fixed_allele, x[0]], in_encl, arg_dict, weights))
	elif ml_params['dimensions'] == 2:
		fn = lambda x: (-weights['frr'] * FRR_sam_likelihood (x, in_frep, arg_dict, weights) + \
						-weights['srp'] * span_sam_likelihood(x, in_span, arg_dict, weights) + \
						-weights['er' ] * encl_sam_likelihood(x, in_encl, arg_dict, weights))

	res = minimize(fn, x0 = ml_params['init_state'], \
						method = 'L-BFGS-B', \
						bounds = ml_params['bounds'])
	return [int(round(j)) for j in res.x]

def likelihood(in_pref, exp_dir, weights, alleles):
	arg_dict = load_profile(exp_dir)
	in_span = in_pref + '_srp.sam'
	in_encl = in_pref + '_er.sam'
	in_frep = in_pref + '_frr.sam'
	fn = lambda x:  -1 * (-weights['frr'] * FRR_sam_likelihood (x, in_frep, arg_dict, weights) + \
					-weights['srp'] * span_sam_likelihood(x, in_span, arg_dict, weights) + \
					-weights['er' ] * encl_sam_likelihood(x, in_encl, arg_dict, weights))
	return fn(alleles)


exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_41_cov60_dist500_DIP/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_45_cov60_dist500_DIP_const40/'

exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_51_cov20_dist500_DIP_const40/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_52_cov60_dist500_DIP_const70/'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_49_cov60_dist500_DIP_const40/'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_48_cov80_dist500_DIP_const40/'
# exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_46_cov80_dist500_DIP_const20/'
# exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_47_cov80_dist500_DIP_const50/'
print exp_dir
weights = {	'frr': 			0.8, \
				'srp': 		1.0, \
				'er': 		1.0, \
				'allele_1': 0.5,\
				'allele_2': 0.5}
oalp_params = {'steps': 	30, \
				'updates': 	0, \
				'Tmax': 	2000.0, \
				'Tmin': 	0.1, \
				'copy_strategy': 'slice', \
				'init_state': [20, 10]}

arg_dict = load_profile(exp_dir)
in_pref_pref = exp_dir + 'aligned_read/nc_'
ml_range = [arg_dict['read_len'] / len(arg_dict['motif']) - 2, 260]
if 'constant_allele' in arg_dict:
	constant_allele = arg_dict['constant_allele']
else:
	constant_allele = arg_dict['ref_allele_count']

print 'const', constant_allele

print '## True\t\tEstimate'
for i in arg_dict['num_copy']:
	# if i > 3 and i <= 30:
	all_alleles = []
	if i > 0:
		in_pref = in_pref_pref + str(i)
		# print in_pref
		# First, genotype enclosing
		allele_list_1, freq_list = encl_sam_genotype(in_pref + '_er.sam', arg_dict)
		all_alleles = all_alleles + allele_list_1
		print allele_list_1, freq_list
		if len(allele_list_1) == 2:
			oalp_params['init_state'] = allele_list_1
		elif len(allele_list_1) == 1:
			oalp_params['init_state'] = allele_list_1 + allele_list_1
		elif len(allele_list_1) > 2:
			oalp_params['init_state'] = top_two(allele_list_1, freq_list)

		# Second, sim annealing for hipstr region
		state, e = run_sim_annealing(exp_dir, in_pref, weights, oalp_params)
		# print state, e
		if state != allele_list_1:
			all_alleles = all_alleles + state

		# Third, run 1-D maximum likelihood for all of current alleles
		for allele in all_alleles:
			ml_params = {'init_state': [ml_range[0] * 2], \
							'bounds': [ml_range], \
							'dimensions': 1, \
							'fixed_allele': allele}
			ml_alleles = run_max_likelihood(in_pref, exp_dir, weights, ml_params)
			# print ml_alleles

			all_alleles = all_alleles + [ml_alleles[0]]

		# Forth, run 2-D maximum likelihood over ml_range
		ml_params = {'init_state': [int(ml_range[0] * 1.5), int(ml_range[0] * 2.5)], \
						'bounds': [ml_range, ml_range], \
						'dimensions': 2}
		ml_alleles = run_max_likelihood(in_pref, exp_dir, weights, ml_params)
		# print ml_alleles

		all_alleles = all_alleles + [ml_alleles[0]]

		print 'All alleles: ', all_alleles
		max_likelihood = -10000
		for a in all_alleles:
			for b in all_alleles:
				like = likelihood(in_pref, exp_dir, weights, [a, b])
				if like > max_likelihood:
					max_likelihood = like
					max_alleles = [a,b]
		print '##', [constant_allele, i], '\t\t', max_alleles


