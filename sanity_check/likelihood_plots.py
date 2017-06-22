import sys, csv,os,errno
import numpy as np
sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile, extract_locus_info
from likelihood_functions import span_genotype_likelihood, FRR_genotype_likelihood, FRR_read_prob, FRR_class_prob
from likelihood_functions import span_allele_likelihood, FRR_allele_likelihood, span_read_prob, span_class_prob, encl_class_prob
from maximum_likelihood_core_function import FRR_sam_likelihood, span_sam_likelihood, ml_FRR_spanning, encl_sam_likelihood
from extract_genome import extract_pre_post_flank

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# StackOverflow
def mkdir_p(path):
    try:
        os.makedirs(path)
        print 'Created path: ', path
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def FRR_read_plot(arg_dict, outPath):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel('Read Probability')
	#ax.set_ylim(left=0, right=100)
	ax.set_xlabel('Omega (Distance from STR region)')
	#ax.set_xlim(bottom=0, top=3000)


	omega_list = range (0, 600, 10)
	allele_len_list = [0, 20, 40, 60, 80, 100, 150, 200, 250]
	colors = cm.rainbow(np.linspace(0, 1, len(allele_len_list)))
	j = 0
	for allele_len in allele_len_list:
		read_prob_omega = []
		for omega in omega_list:
			read_prob_omega.append(FRR_read_prob(arg_dict, allele_len, omega))

		ax.plot(omega_list, read_prob_omega, color = colors[j], lw = 3, label = str(allele_len))
		# ax.plot(xplot, max(y) * gauss(xplot, param1[0], param1[1], param1[2]), color = 'blue', lw = 3, label = 'model1')
		j = j + 1
	ax.legend()
	plt.title('FRR - Read Prob - Mean Ins: '+ str(arg_dict['read_ins_mean']) +' - diff Allele lengths')
	fig.savefig(outPath)

def span_read_plot(arg_dict, outPath):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel('Read Probability')
	#ax.set_ylim(left=0, right=100)
	ax.set_xlabel('Sample Insert Size')
	#ax.set_xlim(bottom=0, top=3000)


	ins_list = range (0, 800, 5)
	allele_len_list = [0, 20, 40, 60, 80, 100, 150, 200, 250]
	colors = cm.rainbow(np.linspace(0, 1, len(allele_len_list)))
	j = 0
	for allele_len in allele_len_list:
		read_prob_ins = []
		for ins in ins_list:
			read_prob_ins.append(span_read_prob(arg_dict, allele_len, ins))

		ax.plot(ins_list, read_prob_ins, color = colors[j], lw = 3, label = str(allele_len))
		# ax.plot(xplot, max(y) * gauss(xplot, param1[0], param1[1], param1[2]), color = 'blue', lw = 3, label = 'model1')
		j = j + 1
	ax.legend()
	plt.title('Spanning - Read Prob - Mean Ins: '+ str(arg_dict['read_ins_mean']) +' - diff Allele lengths')
	fig.savefig(outPath)

def class_plot(arg_dict, sam_pref, outPath_class):
	# allele_len_list = [0,  5,  10,   15,20,   30,40,50,60,70,80,90,100,120,150,180,210,250]
	allele_len_list = arg_dict['num_copy']
	coverage = arg_dict['coverage']
	class_ext = 'er.sam'
	er_class_prob_exp = []
	frr_class_prob_exp = []
	srp_class_prob_exp = []
	er_class_prob_math = []
	frr_class_prob_math = []
	srp_class_prob_math = []
	for allele_len in allele_len_list:
		sam_path_all = sam_pref + str(allele_len) + '.sam'
		count = 0
		total_count = 0.0
		seen_reads = []
		with open(sam_path_all, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						total_count = total_count + 1.0
						seen_reads.append(record[0])

		sam_path_er = sam_pref + str(allele_len) + '_er.sam'
		sam_path_srp = sam_pref + str(allele_len) + '_srp.sam'
		sam_path_frr = sam_pref + str(allele_len) + '_frr.sam'

		er_count = 0.0
		seen_reads = []
		with open(sam_path_er, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						er_count = er_count + 1.0
						seen_reads.append(record[0])

		frr_count = 0.0
		seen_reads = []
		with open(sam_path_frr, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						frr_count = frr_count + 1.0
						seen_reads.append(record[0])

		srp_count = 0.0
		seen_reads = []
		with open(sam_path_srp, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						srp_count = srp_count + 1.0
						seen_reads.append(record[0])

		
		print 'Allele Len:', allele_len
		print 'Total Reads:', total_count
		print 'ER:', er_count, '\tFRR:', frr_count, '\tSRP:', srp_count
		print 'ER:', er_count / total_count, '\tFRR:', frr_count / total_count, '\tSRP:', srp_count / total_count
		print 'ER:', encl_class_prob(arg_dict, allele_len), '\tFRR:', FRR_class_prob(arg_dict, allele_len), '\tSRP:', span_class_prob(arg_dict, allele_len)

		frr_class_prob_exp.append(frr_count / total_count)
		srp_class_prob_exp.append(srp_count / total_count)
		er_class_prob_exp.append(er_count / total_count)

		frr_class_prob_math.append(FRR_class_prob(arg_dict, allele_len))
		srp_class_prob_math.append(span_class_prob(arg_dict, allele_len))
		er_class_prob_math.append(encl_class_prob(arg_dict, allele_len))


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel('Class Probability')
	#ax.set_ylim(left=0, right=100)
	ax.set_xlabel('Allele Length (#copy)')
	#ax.set_xlim(bottom=0, top=3000)

	ax.plot(allele_len_list,frr_class_prob_math, color = 'red', linestyle = '--', lw = 3, label = 'FRR Math')
	ax.plot(allele_len_list,srp_class_prob_math, color = 'blue', linestyle = '--', lw = 3, label = 'Span Math')
	ax.plot(allele_len_list,er_class_prob_math, color = 'green', linestyle = '--', lw = 3, label = 'Enclose Math')

	ax.plot(allele_len_list,frr_class_prob_exp, color = 'red', lw = 3, label = 'FRR Exp')
	ax.plot(allele_len_list,srp_class_prob_exp, color = 'blue', lw = 3, label = 'Span Exp')
	ax.plot(allele_len_list,er_class_prob_exp, color = 'green', lw = 3, label = 'Enclose Exp')

	ax.legend()
	plt.title('Class Prob - Mean Ins: '+ str(arg_dict['read_ins_mean']) +' - Cov' + str(coverage))
	fig.savefig(outPath_class)


def diff_plot(arg_dict, sam_pref, outPath_class):
	# allele_len_list = [0,  5,  10,   15,20,   30,40,50,60,70,80,90,100,120,150,180,210,250]
	allele_len_list = arg_dict['num_copy']
	allele_len_list = [40, 50, 60, 70, 80, 90, 100, 120]
	coverage = arg_dict['coverage']
	class_ext = 'er.sam'
	er_diff_prob_exp = []
	frr_diff_prob_exp = []
	srp_diff_prob_exp = []
	for allele_len in allele_len_list:
		sam_path_all = sam_pref + str(allele_len) + '.sam'
		count = 0
		total_count = 0.0
		seen_reads = []
		with open(sam_path_all, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						total_count = total_count + 1.0
						seen_reads.append(record[0])

		sam_path_er = sam_pref + str(allele_len) + '_er.sam'
		sam_path_srp = sam_pref + str(allele_len) + '_srp.sam'
		sam_path_frr = sam_pref + str(allele_len) + '_frr.sam'

		er_count = 0.0
		seen_reads = []
		with open(sam_path_er, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						er_count = er_count + 1.0
						seen_reads.append(record[0])

		frr_count = 0.0
		seen_reads = []
		with open(sam_path_frr, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						frr_count = frr_count + 1.0
						seen_reads.append(record[0])

		srp_count = 0.0
		seen_reads = []
		with open(sam_path_srp, 'r') as f:
			for record in csv.reader(f, dialect='excel-tab'):
				if record[0][0] != '@':
					if record[0] not in seen_reads:
						srp_count = srp_count + 1.0
						seen_reads.append(record[0])

		
		print 'Allele Len:', allele_len
		print 'Total Reads:', total_count
		print 'ER:', er_count, '\tFRR:', frr_count, '\tSRP:', srp_count
		print 'ER:', er_count / total_count, '\tFRR:', frr_count / total_count, '\tSRP:', srp_count / total_count
		print 'ER:', encl_class_prob(arg_dict, allele_len), '\tFRR:', FRR_class_prob(arg_dict, allele_len), '\tSRP:', span_class_prob(arg_dict, allele_len)

		frr_diff_prob_exp.append(frr_count / total_count - FRR_class_prob(arg_dict, allele_len))
		srp_diff_prob_exp.append(srp_count / total_count - span_class_prob(arg_dict, allele_len))
		er_diff_prob_exp.append(er_count / total_count - encl_class_prob(arg_dict, allele_len))


	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel('Diff Class Probability (real - math)')
	#ax.set_ylim(left=0, right=100)
	ax.set_xlabel('Allele Length (#copy)')
	#ax.set_xlim(bottom=0, top=3000)

	ax.plot(allele_len_list,frr_diff_prob_exp, color = 'red', lw = 3, label = 'FRR Diff')
	ax.plot(allele_len_list,srp_diff_prob_exp, color = 'blue', lw = 3, label = 'Span Diff')
	ax.plot(allele_len_list,er_diff_prob_exp, color = 'green', lw = 3, label = 'Enclose Diff')
	print 'FRR', frr_diff_prob_exp
	print 'SRP', srp_diff_prob_exp
	ax.legend()
	plt.title('Class Prob Diff - Mean Ins: '+ str(arg_dict['read_ins_mean']) +' - Cov' + str(coverage))
	fig.savefig(outPath_class)




def likelihood_plot(arg_dict, sam_pref, weights, outPath):
	coverage = arg_dict['coverage']
	# file_len_list = [0,5,10,15,20,30,40,50,60,70,80,90,100,120,150,180,210,250]
	file_len_list = arg_dict['num_copy']
	file_len_list = [3, 5, 10, 20, 30, 40, 60]
	colors = cm.rainbow(np.linspace(0, 1, len(file_len_list)))

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_ylabel('Likelihood')
	# ax.set_ylim(bottom=-8000, top=0)
	ax.set_xlabel('Allele Length (#copy)')
	#ax.set_xlim(left=0, right=100)
	j = 0
	for file_num in file_len_list:
		in_encl = sam_pref + str(file_num) + '_er.sam'
		in_span = sam_pref + str(file_num) + '_srp.sam'
		in_frep = sam_pref + str(file_num) + '_frr.sam'
		likelihood_array = []
		fixed_allele = 14
		allele_range = range(0, 70)
		for allele in allele_range:
			if 'diploid' in arg_dict and arg_dict['diploid'] == 'True':
				x = [allele, fixed_allele]
			else:
				x = allele
			likelihood = -1 * (-weights['frr'] * FRR_sam_likelihood (x, in_frep, arg_dict, weights) + \
					-weights['srp'] * span_sam_likelihood(x, in_span, arg_dict, weights) + \
					-weights['er' ] * encl_sam_likelihood(x, in_encl, arg_dict, weights))

			likelihood_array.append(likelihood)
		print 'Plotting for file:', file_num
		ax.plot(allele_range,likelihood_array, color = colors[j], lw = 3, label = str(file_num))
		j = j + 1
	ax.legend()
	plt.title('Likelihood - Mean Ins: '+ str(arg_dict['read_ins_mean']) +' - Cov '+str(coverage))
	fig.savefig(outPath)


exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_27_cov60_dist500_hap_viz/'

outPath_FRR_read = '/storage/nmmsv/str-expansions/sanity_check/fig_FRR_read.pdf'
outPath_span_read = '/storage/nmmsv/str-expansions/sanity_check/fig_span_read.pdf'
outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/fig_likelihood.pdf'

sam_pref = '/storage/nmmsv/expansion-experiments/ATXN7_27_cov60_dist500_hap_viz/aligned_read/nc_'
outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_class_60_bug0.pdf'
outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_class_60_true.pdf'
outPath_class_diff = '/storage/nmmsv/str-expansions/sanity_check/diff/ATXN7_60_diff.pdf'

exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_28_cov250_dist500_hap_viz/'
outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_class_250.pdf'
outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_class_250_true.pdf'


exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_30_cov10_dist500_hap_viz/'
outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_10.pdf'
outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_10_true.pdf'

exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/'
outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_60.pdf'
outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_60_true.pdf'
outPath_class_diff = '/storage/nmmsv/str-expansions/sanity_check/diff/ATXN3_60_diff.pdf'

# exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_31_cov10_dist500_hap_viz/'
# outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_CACNA1A_10.pdf'
# outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_CACNA1A_10_true.pdf'


exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_33_cov250_dist500_hap_viz/'
outPath_class = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_250.pdf'
outPath_class_true = '/storage/nmmsv/str-expansions/sanity_check/fig_ATXN3_250_true.pdf'
outPath_class_diff = '/storage/nmmsv/str-expansions/sanity_check/diff/ATXN3_250_diff.pdf'



# sam_pref = exp_dir + 'aligned_read/nc_'
# sam_pref_true = exp_dir + 'aligned_read/true_filter/nc_'
# arg_dict = load_profile(exp_dir)
# class_plot(arg_dict, sam_pref_true, outPath_class_true)
# class_plot(arg_dict, sam_pref, outPath_class)
# diff_plot(arg_dict, sam_pref, outPath_class_diff)


exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_27_cov60_dist500_hap_viz/'
outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/three_class/27_ATXN7_60.pdf'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_28_cov250_dist500_hap_viz/'
# # exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_30_cov10_dist500_hap_viz/'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/'
# # outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/fig_likelihood_wfrr1_0_ATXN3_60.pdf'
# # exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_31_cov10_dist500_hap_viz/'
# # outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/fig_likelihood_CACNA_10.pdf'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_33_cov250_dist500_hap_viz/'
# exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_34_cov250_dist500_hap_viz/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_41_cov60_dist500_DIP/'
outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/three_class/DIP_41_ATXN3_60_95.pdf'
# exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_39_cov60_dist500_hap_viz/'
# outpath_likelihood = '/storage/nmmsv/str-expansions/sanity_check/three_class/39_ATXN3_60.pdf'


arg_dict = load_profile(exp_dir)

# print 'Real\t\tEstimate'
# for i in arg_dict['num_copy']:
# 	file_nc = str(i)
# 	print i,'\t\t', ml_FRR_spanning(exp_dir + 'aligned_read/nc_' + file_nc, exp_dir, weights)
# print 
weights = {	'frr': 0.8, \
			'srp': 1.0, \
			'er': 1.0, \
			'allele_1': 0.5,\
			'allele_2': 0.5}
likelihood_plot(arg_dict, exp_dir + 'aligned_read/nc_', weights, outpath_likelihood)


# outPath_path = exp_dir + 'outfiles/'
# mkdir_p(outPath_path)
# outPath_text_frr_weights = outPath_path + 'frr_weights.txt'
# print outPath_text_frr_weights

# print arg_dict['ref_allele_count']
# handle = open(outPath_text_frr_weights, 'w')
# handle.write('>>>> Trying different weights for each str_len to find the best weight \n')
# for nc in arg_dict['num_copy']:
# 	file_nc = str(nc)
# 	tot_abs_dif = []
# 	print '>> ', nc
# 	handle.write('\n>> ' + file_nc + '\n')
# 	if nc >= 33:
# 		frr_list = [.1,.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2]
# 	else:
# 		frr_list = [1.0]
# 	for frr in frr_list:
# 		weights = {	'frr': frr, \
# 					'srp': 1.0}
# 		print frr,'\t', ml_FRR_spanning(exp_dir + 'aligned_read/nc_' + file_nc, exp_dir, weights)
# 		handle.write(str(frr) + '\t' + str(ml_FRR_spanning(exp_dir + 'aligned_read/nc_' + file_nc, exp_dir, weights)) + '\n')

# handle.close() 




