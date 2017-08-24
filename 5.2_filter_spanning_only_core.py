import subprocess
import os, sys
import numpy as np

sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import expansion_aware_realign, classify_realigned_read, reverse_strand
from load_info import load_profile, extract_locus_info
from extract_genome import extract_pre_post_flank

import argparse

parser = argparse.ArgumentParser('Filter sam/bam files to only keep spanning reads.')
parser.add_argument('--out-pref', 	type = str, required = True)
parser.add_argument('--in-pref', 	type = str, required = True)
parser.add_argument('--exp-dir',	type = str, required = True)

args = parser.parse_args()

out_pref = args.out_pref
in_pref = args.in_pref
exp_dir = args.exp_dir

arg_dict = load_profile(exp_dir)

read_len = arg_dict['read_len']
locus = arg_dict['locus']
motif = arg_dict['motif']

chrom, locus_start, locus_end = extract_locus_info(locus)
pre, post = extract_pre_post_flank(exp_dir, read_len)

score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False
margin = 2


in_sam = in_pref + '.sam'
out_sam = out_pref + '.sam'
out_sam_handle = open(out_sam, 'w')
print 'Filtering ' + in_pref + '.sam'
with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] == '@':
			out_sam_handle.write(record)
		else:
			row = record.split()
			# According to SAM specification:
			RNAME = row[2]
			POS = int(row[3])
			RNEXT = row[6]
			PNEXT = int(row[7])
			TLEN = int(row[8])
			SEQ = row[9]

			if (RNAME == chrom and POS <= locus_start - read_len and RNEXT == '=' and PNEXT >= locus_end):
				print TLEN
				nc_col = 'nc:i:' + str(0)
				ps_col = 'ps:i:' + str(0)
				sc_col = 'sc:i:' + str(0)
				rc_col = 'rc:Z:' + 'unknown'
				is_col = 'is:i:' + str(int(np.abs(TLEN)))
				out_sam_handle.write('\t'.join(row + [nc_col, ps_col, sc_col, rc_col, is_col])+'\n')
			elif (RNAME == chrom and POS <= locus_start - read_len and RNEXT == '=' and PNEXT <= locus_start) or \
				(RNAME == chrom and POS >= locus_end and RNEXT == '=' and PNEXT >= locus_end):
				# ignoring read pairs that are both before or after STR region.
				pass
			else:
				# We realign first read to check if it's pre or post flanking
				nCopy, pos, score = expansion_aware_realign(SEQ, pre, post, motif, score_dict, verbose)
				nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(SEQ), pre, post, motif, score_dict, verbose)
				if score_rev > score:
					nCopy = nCopy_rev
					pos = pos_rev
					score = score_rev
					SEQ = reverse_strand(SEQ)
				read_class = classify_realigned_read(SEQ, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
				# We're allowing flanking from both sides, but one side should be mapped to the correct location for this to work.
				# i.e. we don't rescue the case that both reads are flanking, and both are mapped to a different locus.
				nc_col = 'nc:i:' + str(nCopy)
				ps_col = 'ps:i:' + str(pos)
				sc_col = 'sc:i:' + str(score)
				rc_col = 'rc:Z:' + read_class

				if ((RNAME == chrom and RNEXT == '=') or RNEXT == chrom) and PNEXT >= locus_end - read_len and read_class =='PreFlank':
					abs_insert_size = np.abs((locus_start + nCopy * len(motif) - read_len) - PNEXT) + read_len
					print 'Found ', read_class, '\t', abs_insert_size, '\t', pos, '\t', row[0]	
					is_col = 'is:i:' + str(int(abs_insert_size) - 1) # -1 to match CPP
					out_sam_handle.write('\t'.join(row + [nc_col, ps_col, sc_col, rc_col, is_col])+'\n')
					print
				elif ((RNAME == chrom and RNEXT == '=') or RNEXT == chrom) and PNEXT <= locus_start and read_class =='PostFlank':
					
					abs_insert_size = np.abs((locus_end - nCopy * len(motif)) - PNEXT) + read_len
					print 'Found ', read_class, '\t', abs_insert_size, '\t', pos, '\t', row[0]
					is_col = 'is:i:' + str(int(abs_insert_size) + 1) # +1 to match CPP
					out_sam_handle.write('\t'.join(row + [nc_col, ps_col, sc_col, rc_col, is_col])+'\n')
					print

out_sam_handle.close()

# Sort and index filtered bam
# os.system('samtools view -bT ' + \
# 		ref_gen_dir + ' ' + \
# 		out_pref + '.sam' + ' ' + 
# 		'> ' + out_pref + '.bam')
# os.system('samtools sort -o ' + \
# 		out_pref + '.sorted.bam '+ \
# 		out_pref + '.bam')
# os.system('samtools index ' + \
# 		out_pref + '.sorted.bam ' + \
# 		out_pref + '.sorted.bai')
# os.system('samtools index ' + out_pref + '.sorted.bam')
