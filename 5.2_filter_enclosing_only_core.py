import subprocess
import os, sys
import argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import expansion_aware_realign, classify_realigned_read, reverse_strand
from load_info import load_profile, extract_locus_info
from extract_genome import extract_pre_post_flank

def realignment_string(temp_fa_dir, read_len, motif):
	with open(temp_fa_dir + '_prefix.fa', 'r') as f:
		f.readline()
		pref = f.readline()[-read_len - 1:].strip()
	with open(temp_fa_dir + '_suffix.fa', 'r') as f:
		f.readline()
		suff = f.readline()[:read_len].strip()
	str_region = (int(read_len / len(motif)) + 1) * motif
	# str_region = str_region [0 : read_len]
	return pref + str_region + suff



parser = argparse.ArgumentParser('Filter sam/bam files to only keep reads with one repeating mate.')
parser.add_argument('--out-pref', 	type = str, required = True)
parser.add_argument('--in-pref', 	type = str, required = True)
parser.add_argument('--exp-dir',	type = str, required = True)


args = parser.parse_args()

out_pref = args.out_pref
in_pref = args.in_pref
exp_dir = args.exp_dir

arg_dict = load_profile(exp_dir)

motif = arg_dict['motif']
read_len = arg_dict['read_len']
locus = arg_dict['locus']
ref_gen_dir = arg_dict['ref_gen_dir']
exp_name = arg_dict['exp_name']
temp_dir = exp_dir + '/temp/'
temp_fa_dir = temp_dir + exp_name
read_ins_mean = arg_dict['read_ins_mean']
read_ins_stddev = arg_dict['read_ins_stddev']

chrom, locus_start, locus_end = extract_locus_info(locus)
pre, post = extract_pre_post_flank(exp_dir, read_len)

score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False
margin = 2

vicinity_factor = 3
far_factor = 6

# Setting the filter edges for left and right
# Looking for reads such that r1 <= left_filter AND r2 >= right_filter
left_filter = locus_start - read_len
right_filter = locus_end

realign = realignment_string (temp_fa_dir, read_len, motif)

in_sam = in_pref + '.sam'
print 'Filtering ' + in_pref + '.sam'

list_reads = []
out_sam = out_pref + '.sam'
out_sam_handle = open(out_sam, 'w')
with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] != '@':
			row = record.split()
			# According to SAM specification:
			RNAME = row[2]
			POS = int(row[3])
			RNEXT = row[6]
			PNEXT = int(row[7])
			TLEN = int(row[8])
			SEQ = row[9]

			if len(SEQ) > 0.9 * read_len:		# Discarding soft clipped reads
				# 1) if read is overlapping with STR region
				if (RNAME == chrom and POS >= locus_start - 2 * read_len and POS <= locus_start + read_len):
					# Performing realignment to check if IRR
					nCopy, pos, score = expansion_aware_realign(SEQ, pre, post, motif, score_dict, verbose)
					nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(SEQ), pre, post, motif, score_dict, verbose)
					if score_rev > score:
						nCopy = nCopy_rev
						pos = pos_rev
						score = score_rev
						SEQ = reverse_strand(SEQ)
					read_class = classify_realigned_read(SEQ, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
					if read_class == 'Enclosing' or\
						 (nCopy == 0 and score > 0.95 * read_len * score_dict['match']and POS > 0.2 * read_len and POS < 0.8 * read_len):
						nc_col = 'nc:i:' + str(nCopy)
						ps_col = 'ps:i:' + str(pos)
						sc_col = 'sc:i:' + str(score)
						rc_col = 'rc:Z:' + read_class
						out_sam_handle.write('\t'.join(row + [nc_col, ps_col, sc_col, rc_col])+'\n')
						print 'Found Enclosing Read!'
						print nCopy, pos, score, SEQ
					else:
						pass
				# 2) Mate/next is in vicinity of STR region, but read is not aligned to the same chrom (ins = 0)
				# 			OR aligned way too far (abs(ins) > read_ins_mean + 6 * read_ins_stddev)
				elif (RNEXT == chrom and \
						PNEXT >= locus_start - read_ins_mean - vicinity_factor * read_ins_stddev and \
						PNEXT <= locus_end   + read_ins_mean + vicinity_factor * read_ins_stddev) and \
						(TLEN == 0 or np.abs(TLEN) > read_ins_mean + far_factor * read_ins_stddev):
					# Performing realignment to check if IRR
					nCopy, pos, score = expansion_aware_realign(SEQ, pre, post, motif, score_dict, verbose)
					nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(SEQ), pre, post, motif, score_dict, verbose)
					if score_rev > score:
						nCopy = nCopy_rev
						pos = pos_rev
						score = score_rev
						SEQ = reverse_strand(SEQ)
					read_class = classify_realigned_read(SEQ, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
					if read_class == 'Enclosing' or \
						(nCopy == 0 and score > 0.95 * read_len * score_dict['match'] and POS > 0.2 * read_len and POS < 0.8 * read_len):
						nc_col = 'nc:i:' + str(nCopy)
						ps_col = 'ps:i:' + str(pos)
						sc_col = 'sc:i:' + str(score)
						rc_col = 'rc:Z:' + read_class
						out_sam_handle.write('\t'.join(row + [nc_col, ps_col, sc_col, rc_col])+'\n')
						print 'Found Enclosing Read!'
						print nCopy, pos, score, SEQ
					else:
						pass
		else:
			out_sam_handle.write(record)


out_sam_handle.close()

# print list_rows
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
