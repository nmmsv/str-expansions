from likelihood_spanning import calc_likelihood_span
import csv, sys
import numpy as np
from scipy.optimize import minimize
import pickle, argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile, extract_locus_info
from likelihood_functions import span_genotype_likelihood, IRR_genotype_likelihood

def IRR_sam_likelihood (A, B, sam_path, arg_dict):
	locus = arg_dict['locus']
	chrom, locus_start, locus_end = extract_locus_info(locus)
	log_likelihood = 1
	with open(sam_path, 'r') as irr_handle:
		for record in csv.reader(irr_handle, dialect = 'excel-tab'):
			if record[0][0] != '@':
				sample_dfl = locus_start - int(record[3])
				samp_likelihood = IRR_genotype_likelihood(arg_dict, A, B, sample_dfl)
				# print record[0], '\t', locus_start - int(record[3]), '\t', samp_likelihood
				if samp_likelihood > 0:
					samp_log_likelihood = np.log(samp_likelihood)
				else:
					samp_log_likelihood = -500
				log_likelihood = log_likelihood + samp_log_likelihood
				# log_likelihood = log_likelihood * samp_likelihood
	return log_likelihood


def span_sam_likelihood (A, B, sam_path, arg_dict):
	locus = arg_dict['locus']
	chrom, locus_start, locus_end = extract_locus_info(locus)
	log_likelihood = 0
	with open(sam_path, 'r') as irr_handle:
		for record in csv.reader(irr_handle, dialect = 'excel-tab'):
			if record[0][0] != '@':
				sample_ins = int(record[8])
				samp_likelihood = span_genotype_likelihood(arg_dict, A, B, sample_ins)

				print sample_ins, samp_likelihood
				# print record[0], '\t', int(record[8]), '\t', samp_likelihood
				if samp_likelihood > 0:
					samp_log_likelihood = np.log(samp_likelihood)
				else:
					samp_log_likelihood = -500
				log_likelihood = log_likelihood + samp_log_likelihood
	return log_likelihood



parser = argparse.ArgumentParser('Use filtered sam files to compute maximum likelihood genotype')
parser.add_argument('--in-sam-span', 	type = str, required = True)
parser.add_argument('--in-sam-IRR', 	type = str, required = True)
parser.add_argument('--exp-dir',	type = str, required = True)


args = parser.parse_args()

span_path = args.in_sam_span
IRR_path  = args.in_sam_IRR
exp_dir = args.exp_dir

arg_dict = load_profile(exp_dir)


# fn = lambda x: -1 * calc_likelihood_sam(x[0], x[1], sam_path, arg_dict)

A = 10
B = 400



fn = lambda x: (-1 * IRR_sam_likelihood(x[0], x[1], IRR_path, arg_dict) + \
				-1 * span_sam_likelihood(x[0], x[1], span_path, arg_dict))

fn_span = lambda x: (-1 * span_sam_likelihood(x[0], x[1], span_path, arg_dict))

fn_IRR = lambda x: (-1 * IRR_sam_likelihood(x[0], x[1], IRR_path, arg_dict))


print fn_span([200, 210])

# # Try dioffenet methods to see which one works best/fastest
# res = minimize(fn, x0 = [10, 40], \
# 				method = 'L-BFGS-B', \
# 				bounds = ((1, 250), (1, 250)))	
# print [round(j) for j in res.x]
# print fn(res.x)
# print fn([10, 80])


res = minimize(fn_IRR, x0 = [210,210], \
				method = 'L-BFGS-B', \
				bounds = ((1, 250), (1, 250)))	
print [round(j) for j in res.x]
print fn_IRR(res.x)
print fn_IRR([250, 250])

# with open ('/storage/nmmsv/expansion-experiments/span.csv', 'w') as csvSpan:
# 	with open ('/storage/nmmsv/expansion-experiments/IRR.csv', 'w') as csvIRR:
# 		with open ('/storage/nmmsv/expansion-experiments/both.csv', 'w') as csvBoth:
# 			writerSpan = csv.writer(csvSpan, delimiter = ',')
# 			writerIRR = csv.writer(csvIRR, delimiter = ',')
# 			writerBoth = csv.writer(csvBoth, delimiter = ',')
# 			for i in range(1, 250):
# 				lstSpan = []
# 				lstIRR = []
# 				lstBoth = []
# 				for j in range(1, 250):
# 					sys.stdout.write("Progress: %d, %d   \r" % (i, j) )
# 					sys.stdout.flush()
# 					spn = fn_span([i, j])
# 					irr = fn_IRR([i, j])
# 					bth = spn + irr
# 					lstSpan.append(spn)
# 					lstIRR.append(irr)
# 					lstBoth.append(bth)
# 				print
# 				print 'Row added!'
# 				writerSpan.writerow(lstSpan)
# 				writerIRR.writerow(lstIRR)
# 				writerBoth.writerow(lstBoth)

# with open ('/storage/nmmsv/expansion-experiments/IRR.csv', 'w') as csvfile:
# 	filewriter = csv.writer(csvfile, delimiter = ',')
# 	for i in range(1, 250):
# 		lst = []
# 		for j in range(1, 250):
# 			lst.append(fn_IRR([i, j]))
# 		print 'Row added! (IRR)'
# 		filewriter.writerow(lst)

# with open ('/storage/nmmsv/expansion-experiments/both.csv', 'w') as csvfile:
# 	filewriter = csv.writer(csvfile, delimiter = ',')
# 	for i in range(1, 250):
# 		lst = []
# 		for j in range(1, 250):
# 			lst.append(fn([i, j]))
# 		print 'Row added! (both)'
# 		filewriter.writerow(lst)

# for i in range(10, 210, 20):
# 	print fn_span([10, i])


