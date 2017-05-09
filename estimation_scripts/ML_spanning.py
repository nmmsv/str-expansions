from likelihood_spanning import calc_likelihood_span
import csv
import numpy as np
from scipy.optimize import minimize
import pickle
# TODO:
# you still need to pass some parameters (not wgsim parameters)
# wgsim parameters will be retrieved from the profile, so pass the path
# to profile as well.
real_alt_allele = 80
exp_name = 'ATXN7_16_testML_cov40_dist1000'
exp_base = '/storage/nmmsv/expansion-experiments/'
exp_dir = exp_base + exp_name

sam_path = exp_dir + '/aligned_read/nc_' + str(real_alt_allele) + '_flt.sam'

try:
	with open(exp_dir + '/profile.txt', 'r') as f:
		arg_dict = pickle.load(f)
except:
	raise
print arg_dict

def calc_likelihood_sam (A, B, sam_path, arg_dict):
	log_like = 0
	with open(sam_path, 'r') as sam_file_handle:
		j = 0
		for line in csv.reader(sam_file_handle, dialect = 'excel-tab'):
			if line[0][0] != '@':
				samp_read = int(line[8])
				log_like = log_like + np.log(calc_likelihood_span(arg_dict, A, B, samp_read))
	return log_like


fn = lambda x: -1 * calc_likelihood_sam(x[0], x[1], sam_path, arg_dict)

# # Try dioffenet methods to see which one works best/fastest
# res = minimize(fn, x0 = [10, 20], \
# 				method = 'L-BFGS-B', \
# 				bounds = ((1, 150), (1, 150)))	

# print [round(j) for j in res.x]
# print fn(res.x)

# for i in range(10, 210, 20):
# 	print fn([10, i])

print fn([10, 20])