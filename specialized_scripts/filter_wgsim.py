import sys, os, errno
sys.path.append('/storage/nmmsv/str-expansions/functions/')
from load_info import load_profile, extract_locus_info

# import argparse
# parser = argparse.ArgumentParser('Filter sam solely based on wgsim coordinates.')
# parser.add_argument('--align', type = str, required = True)
# args = parser.parse_args()
# align_flag = args.align

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



exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_27_cov60_dist500_hap_viz/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN7_28_cov250_dist500_hap_viz/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_30_cov10_dist500_hap_viz/'
exp_dir = '/storage/nmmsv/expansion-experiments/ATXN3_32_cov60_dist500_hap_viz/'
# exp_dir = '/storage/nmmsv/expansion-experiments/CACNA1A_31_cov10_dist500_hap_viz/'



sam_dir = exp_dir + 'aligned_read/'
true_filt_dir = exp_dir + 'aligned_read/true_filter/'
# true_filt_dir = '/storage/nmmsv/python_playground/test_filter_IRR/'
print 'kek'
print true_filt_dir
mkdir_p(true_filt_dir)

arg_dict = load_profile(exp_dir)
flank_len = arg_dict['flank_len']
read_len = arg_dict['read_len']
locus = arg_dict['locus']
motif = arg_dict['motif']
chrom, locus_start_ref, locus_end_ref = extract_locus_info(locus)



nc_list = [0,5,10,15,20,30,40,50,60,70,80,90,100,120,150,180,210,250]

nc_list = [0,3,5,7,10,12,15,20,25,30,40,50,60,70,80,90,100,120,150,180,210,250]

for nCopy in nc_list:
	locus_start = locus_start_ref
	locus_end = locus_start_ref + nCopy * len(motif)

	frr_handle = open(true_filt_dir + 'nc_' + str(nCopy) + '_frr.sam', 'w')
	frr_bank = []
	er_handle = open(true_filt_dir + 'nc_' + str(nCopy) + '_er.sam', 'w')
	srp_handle = open(true_filt_dir + 'nc_' + str(nCopy) + '_srp.sam', 'w')
	sam_copy_handle = open(true_filt_dir + 'nc_' + str(nCopy) + '.sam', 'w')
	sam_file = sam_dir + 'nc_' + str(nCopy) + '.sam'
	with open(sam_file, 'r') as sam_handle:
		for record in sam_handle:
			sam_copy_handle.write(record)
			row = record.split()

			if row[0][0] != '@':
				# According to SAM specification:
				QNAME = row[0]
				RNAME = row[2]
				TLEN = int(row[8])
				SEQ = row[9]

				if len(SEQ) >= 0.9 * read_len:
					QNAME_break = QNAME.split('_')
					Q1 = locus_start_ref - flank_len + int(QNAME_break[8])
					Q2 = locus_start_ref - flank_len + int(QNAME_break[9]) - read_len + 1
					if (Q1 < locus_start and Q1 + read_len > locus_end and TLEN > 0) or \
							(Q2 < locus_start and Q2 + read_len > locus_end and TLEN < 0):
						er_handle.write(record)
					if (Q1 >= locus_start and Q1 + read_len <= locus_end) or \
							(Q2 >= locus_start and Q2 + read_len <= locus_end) and RNAME == chrom:
						frr_handle.write(record)
					if (Q1 <= locus_start and Q2 >= max(locus_end - read_len, locus_start)) or \
							(Q2 <= locus_start and Q1 >= max(locus_end - read_len, locus_start)):
						srp_handle.write(record)
					# if Q1 > locus_start_ref - read_len and Q1 
					# print Q1, Q2
					# print row[3], row[7]
					# print 
			else:
				er_handle.write(record)
				frr_handle.write(record)
				srp_handle.write(record)