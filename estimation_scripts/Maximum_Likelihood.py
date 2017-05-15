from likelihood_spanning import calc_likelihood_span
import csv, sys
import numpy as np
from scipy.optimize import minimize
import pickle, argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile, extract_locus_info
from likelihood_functions import span_genotype_likelihood, IRR_genotype_likelihood
from realignment import expansion_aware_realign
from extract_genome import extract_pre_post_flank
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
					samp_log_likelihood = -250
				log_likelihood = log_likelihood + samp_log_likelihood
				# log_likelihood = log_likelihood * samp_likelihood
	return log_likelihood


def span_sam_likelihood (A, B, sam_path, arg_dict):
	locus = arg_dict['locus']
	chrom, locus_start, locus_end = extract_locus_info(locus)
	log_likelihood = 0
	with open(sam_path, 'r') as irr_handle:
		for record in csv.reader(irr_handle, dialect = 'excel-tab'):
			if record[0][0] != '@' and int(record[8]) != 0:
				sample_ins = int(record[8])
				samp_likelihood = span_genotype_likelihood(arg_dict, A, B, sample_ins)

				# print record[0], '\t', int(record[8]), '\t', np.abs(int(record[8]) - 500) / 3 + 10, '\t', samp_likelihood
				if samp_likelihood > 0:
					samp_log_likelihood = np.log(samp_likelihood)
				elif np.abs(samp_likelihood) < 10**-20:		# accounting for comutational errors
					samp_log_likelihood = np.log(np.abs(samp_likelihood))
				elif samp_likelihood == 0:
					samp_log_likelihood = -250
				else:
					print 'Error! Negative likelihood:', samp_likelihood
				log_likelihood = log_likelihood + samp_log_likelihood
	return log_likelihood

def encl_sam_genotype (sam_path, arg_dict):
	exp_dir = arg_dict['exp_dir']
	read_len = arg_dict['read_len']
	motif = arg_dict['motif']
	score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
	verbose = False	
	pre, post = extract_pre_post_flank(exp_dir, read_len)

	nCopy_dict = {}
	total_count = 0
	with open(sam_path, 'r') as encl_handle:
		for record in csv.reader(encl_handle, dialect = 'excel-tab'):
			if record[0][0] != '@':
				sample = record[9]
				nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
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



parser = argparse.ArgumentParser('Use filtered sam files to compute maximum likelihood genotype')
# parser.add_argument('--in-sam-span', 	type = str, required = True)
# parser.add_argument('--in-sam-IRR', 	type = str, required = True)
parser.add_argument('--in-pref',		type = str, required = True)
parser.add_argument('--exp-dir',	type = str, required = True)



args = parser.parse_args()

# span_path = args.in_sam_span
# IRR_path  = args.in_sam_IRR
in_pref = args.in_pref
exp_dir = args.exp_dir

arg_dict = load_profile(exp_dir)

in_span = in_pref + '_srp.sam'
in_encl = in_pref + '_er.sam'
in_frep = in_pref + '_frr.sam'

# fn = lambda x: -1 * calc_likelihood_sam(x[0], x[1], sam_path, arg_dict)

A = 10
B = 400



fn = lambda x: (-1 * IRR_sam_likelihood(x[0], x[1], in_frep, arg_dict) + \
				-1 * span_sam_likelihood(x[0], x[1], in_span, arg_dict))

fn_span = lambda x: (-1 * span_sam_likelihood(x[0], x[1], in_span, arg_dict))

fn_frep = lambda x: (-1 * IRR_sam_likelihood(x[0], x[1], in_frep, arg_dict))





# encl_allele,encl_allele_freq = encl_sam_genotype(in_encl, arg_dict)
# print encl_allele, encl_allele_freq


encl_allele = [10]
if len(encl_allele) == 1:
	fn_span_single = lambda x: (-1 * span_sam_likelihood(x[0], encl_allele[0], in_span, arg_dict))

	# def diminisher(x, encl_allele):
	# 	val = fn_span_single(encl_allele)
	# 	sig = 1
	# 	amplitude_damper = 0.01
	# 	return amplitude_damper * val * np.exp(-0.5 * ((x - encl_allele[0]) / sig) ** 2)
	# fn_span_single_diminished = lambda x: fn_span_single(x) + diminisher(x[0], encl_allele)
	# Now try using spanning read pairs:
	res = minimize(fn_span_single, x0 = [30], \
					method = 'L-BFGS-B', \
					bounds = [(20, 250)])
	print [round(j) for j in res.x]
	print fn_span_single(res.x)

	
	# with open ('/storage/nmmsv/expansion-experiments/fn_span_single.txt', 'w') as f:
	# 	for i in range(1,100,5):
	# 		func_val = fn_span_single([i])
	# 		print i, func_val
	# 		f.write(str(i) + '\t' + str(func_val) + '\n')



# # Try dioffenet methods to see which one works best/fastest
# res = minimize(fn, x0 = [10, 40], \
# 				method = 'L-BFGS-B', \
# 				bounds = ((1, 250), (1, 250)))	
# print [round(j) for j in res.x]
# print fn(res.x)
# print fn([10, 80])


# res = minimize(fn_IRR, x0 = [210,210], \
# 				method = 'L-BFGS-B', \
# 				bounds = ((1, 250), (1, 250)))	
# print [round(j) for j in res.x]
# print fn_IRR(res.x)
# print fn_IRR([250, 250])

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


