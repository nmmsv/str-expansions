import sys
sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import expansion_aware_realign, classify_realigned_read
from load_info import load_profile, extract_locus_info
from extract_genome import extract_pre_post_flank

read_class = 'srp'
nCopy = 70
filt_path = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/aligned_read/nc_'+str(nCopy)+'_'+read_class+'.sam'
# filt_path = '/storage/nmmsv/python_playground/test_filter_IRR/nc_'+str(nCopy)+'.sam'
filt_path_true = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/aligned_read/true_filter/nc_'+str(nCopy)+'_'+read_class+'.sam'
sam_path = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/aligned_read/nc_'+str(nCopy)+'.sam'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/'
arg_dict = load_profile(exp_dir)
locus = arg_dict['locus']
read_len = arg_dict['read_len']
motif = arg_dict['motif']
chrom, locus_start_ref, locus_end_ref = extract_locus_info(locus)
pre, post = extract_pre_post_flank(exp_dir, read_len)
score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False
margin = 2

print locus_start_ref, locus_end_ref
true_reads = []
kk = 0
with open (filt_path, 'r') as handle:
	for record in handle:
		if record[0] != '@':
			kk = kk + 1
			QNAME = record.split()[0]
			true_reads.append(QNAME)
ll = 0
with open (filt_path_true, 'r') as handle:
	for record in handle:
		if record[0] != '@':
			QNAME = record.split()[0]
			SEQ = record.split()[9]
			ll = ll + 1
			if QNAME not in true_reads:
				if QNAME == 'ATXN7_27_cov60_dist500_hap_viz_50_haplo_2617_3124_0:0:0_0:0:0_f':
					print
					print
				print record
print kk, ll
				# nCopy, pos, score = expansion_aware_realign(SEQ, pre, post, motif, score_dict, verbose)
				# read_class = classify_realigned_read(SEQ, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
				# if read_class == 'IRR':
				# 	print nCopy, score
				# 	print record
				# 	print 
