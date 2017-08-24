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


vicinity_factor = 4

in_sam = in_pref + '.sam'
print 'Filtering ' + in_pref + '.sam'

list_reads = []
out_sam_new = out_pref + '_new.sam'
out_sam_handle_new = open(out_sam_new, 'w')
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
				if (RNAME == chrom and RNEXT == '=') or (RNEXT == chrom):			# checking correct chrom for mate 
					# checking if mate is in vicinity of STR region
					if (PNEXT <= locus_start - read_len and PNEXT >= locus_start - read_len - read_ins_mean - vicinity_factor * read_ins_stddev) or \
						(PNEXT >= locus_end and PNEXT <= locus_end + read_ins_mean + vicinity_factor * read_ins_stddev):
						if row[2] != chrom:
							# Performing realignment to check if IRR
							sample = row[9]
							nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
							nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(sample), pre, post, motif, score_dict, verbose)
							if score_rev > score:
								nCopy = nCopy_rev
								pos = pos_rev
								score = score_rev
								sample = reverse_strand(sample)
							read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
							if read_class == 'IRR':
								list_reads.append(row[0])			# Saving read_ID of potential reads that going to be extracted
								if PNEXT <= locus_start:
									om_col = 'om:i:' + str(locus_start - PNEXT - read_len + 1)
								elif PNEXT >= locus_end:
									om_col = 'om:i:' + str(PNEXT - locus_end - 1)
								else:
									print 'Weird PNEXT:', PNEXT
								out_sam_handle_new.write('\t'.join(row + [om_col])+'\n')
								print 'Found IRR!'
							else:
								pass
						elif int(row[3]) >= locus_start and int(row[3]) <= locus_end - read_len:
							# Performing realignment to check if IRR
							sample = row[9]
							nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
							nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(sample), pre, post, motif, score_dict, verbose)
							if score_rev > score:
								nCopy = nCopy_rev
								pos = pos_rev
								score = score_rev
								sample = reverse_strand(sample)
							read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
							if read_class == 'IRR':
								list_reads.append(row[0])			# Saving read_ID of potential reads that going to be extracted
								if PNEXT <= locus_start:
									om_col = 'om:i:' + str(locus_start - PNEXT - read_len + 1)
								elif PNEXT >= locus_end:
									om_col = 'om:i:' + str(PNEXT - locus_end - 1)
								else:
									print 'Weird PNEXT:', PNEXT
								out_sam_handle_new.write('\t'.join(row + [om_col])+'\n')
							else:
								pass

out_sam = out_pref + '.sam'
out_sam_handle = open(out_sam, 'w')
with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] == '@':
			out_sam_handle.write(record)
		else:
			row = record.split()
			RNAME = row[2]
			POS = int(row[3])
			RNEXT = row[6]
			PNEXT = int(row[7])
			TLEN = int(row[8])
			SEQ = row[9]
			if len(SEQ) > 0.9 * read_len:		# Discarding soft clipped reads
				if row[0] in list_reads:
					# perform realignment to check this read is the non-spanning read (not the IRR mate)
					nCopy, pos, score = expansion_aware_realign(SEQ, pre, post, motif, score_dict, verbose)
					nCopy_rev, pos_rev, score_rev = expansion_aware_realign(reverse_strand(SEQ), pre, post, motif, score_dict, verbose)
					if score_rev > score:
						nCopy = nCopy_rev
						pos = pos_rev
						score = score_rev
						SEQ = reverse_strand(SEQ)
					read_class = classify_realigned_read(SEQ, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
					# print read_class, ':  ', RNAME
					# print sample
					print 
					if read_class == 'NoSpan':
						if POS <= locus_start:
							om_col = 'om:i:' + str(locus_start - POS - read_len + 1)
						elif POS >= locus_end:
							om_col = 'om:i:' + str(POS - locus_end - 1)
						else:
							print 'Weird POS:', POS
						if row[0] == 'ATXN7_27_cov60_dist500_hap_viz_50_haplo_2617_3124_0:0:0_0:0:0_f':
							print 'Gotchaaaa'
						out_sam_handle.write('\t'.join(row + [om_col])+'\n')
						print 'Record Added!'



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
