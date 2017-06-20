import sys, csv,os,errno
import numpy as np
sys.path.append('/storage/nmmsv/str-expansions/functions/')
sys.path.append('/storage/nmmsv/str-expansions/estimation_scripts/')
from load_info import load_profile, extract_locus_info
from likelihood_functions import span_genotype_likelihood, FRR_genotype_likelihood, FRR_read_prob, FRR_class_prob
from likelihood_functions import span_allele_likelihood, FRR_allele_likelihood, span_read_prob, span_class_prob, encl_class_prob
from maximum_likelihood_core_function import FRR_sam_likelihood, span_sam_likelihood, ml_FRR_spanning, ml_encl_FRR_span
from extract_genome import extract_pre_post_flank


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


main_dir = '/storage/nmmsv/expansion-experiments/'


exp_dir = main_dir + 'rpt_ATXN3_1_nC50_cov60_dist500/'
outPath_path = exp_dir + 'outfiles/'
mkdir_p(outPath_path)


outPath_text_frr_weights = outPath_path + 'repeated_experiment.txt'


# print 'Real\t\tEstimate'
# for i in arg_dict['num_copy']:
# 	file_nc = str(i)
# 	print i,'\t\t', ml_FRR_spanning(exp_dir + 'aligned_read/nc_' + file_nc, exp_dir, weights)
# print 

# likelihood_plot(arg_dict, exp_dir + 'aligned_read/nc_', weights, outpath_likelihood)

arg_dict = load_profile(exp_dir)
handle = open(outPath_text_frr_weights, 'w')
frr_weight = 0.8
handle.write('>>>> Setting FRR weight to '+str(frr_weight)+' and finding the max likelihood estimate for repeated experiments.\n')
weights = {	'frr': frr_weight, \
			'srp': 1.0, \
			'er' : 1.0}	
j = 0
print '>> ', arg_dict['num_copy']
handle.write('\n>> ' + str(arg_dict['num_copy']) + '\n')
result_array = []
for nc in arg_dict['num_copy']:
	file_nc = str(nc)
	# result = ml_FRR_spanning(exp_dir + 'aligned_read/'+str(j)+'nc_' + file_nc, exp_dir, weights)
	result = ml_encl_FRR_span(exp_dir + 'aligned_read/'+str(j)+'nc_' + file_nc, exp_dir, weights)
	print j,'\t', result[0]
	handle.write(str(j) + '\t' + str(result[0]) + '\n')
	result_array.append(result[0])
	j = j + 1

print '--------'
print 'Estimate Mean: ' +  str(np.mean(result_array)) + '\n'
print 'Estimate Standard Deviation: ' + str(np.std(result_array)) + '\n'
handle.write('\n-------------\n')
handle.write('Estimate Mean: ' +  str(np.mean(result_array)) + '\n')
handle.write('Estimate Standard Deviation: ' + str(np.std(result_array)) + '\n')
handle.close() 
