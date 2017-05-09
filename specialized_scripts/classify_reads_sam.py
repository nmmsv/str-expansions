import argparse, sys

sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import expansion_aware_realign, classify_realigned_read
from load_info import load_profile, extract_locus_info
from extract_genome import extract_pre_post_flank


parser = argparse.ArgumentParser('A test script to classify all reads in a sam file.')
parser.add_argument('--in-sam', 	type = str, required = True)
parser.add_argument('--exp-dir',	type = str, required = True)

args = parser.parse_args()

in_sam = args.in_sam
exp_dir = args.exp_dir


arg_dict = load_profile(exp_dir)

read_len = arg_dict['read_len']
motif = arg_dict['motif']
locus_bed = arg_dict['locus']

pre, post = extract_pre_post_flank(exp_dir, read_len)
score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False
margin = 2
chrom, locus_start, locus_end = extract_locus_info(locus_bed)

with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] == '@':
			# out_sam_handle.write(record)
			pass
		else:
			row = record.split()

			if (row[2] == chrom and int(row[3]) >= locus_start - read_len and int(row[3]) <= locus_end): # Reads mapped within region
				sample = row[9]
				nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
				read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
				print '>> ', read_class, ':'
				print '>> nCopy =', nCopy
				print '>> score =', score
				print '>> pos =', pos
				print sample
				print

			# if (row[2] == chrom and row[6] == '=') or (row[6] == chrom):		# checking correct chrom for mate 
			# 	if int(row[7]) <= locus_start:									# checking if mate is before STR region
			# 		if row[2] != chrom:
			# 			out_sam_handle.write(record)
			# 			print row[9]
			# 			pos = smith_waterman(realign, row[9], score_dict, verbose)
			# 			if pos > read_len - 2 and pos < read_len + 4:
			# 				print '## IRR'
			# 			elif pos <= read_len - 2:
			# 				print '## Pre-Flank', row[0]
			# 			else:
			# 				print '## Post-Flank (mate before STR)', row[0]
			# 		elif int(row[3]) >= locus_start - read_len and int(row[3]) <= locus_end:
			# 			# list_rows.append(row[0])

			# 			out_sam_handle.write(record)
			# 			print row[9]
			# 			pos = smith_waterman(realign, row[9], score_dict, verbose)
			# 			if pos > read_len - 2 and pos < read_len + 4:
			# 				print '## IRR'
			# 			elif pos <= read_len - 2:
			# 				print '## Pre-Flank', row[0]
			# 			else:
			# 				print '## Post-Flank (mate before STR)', row[0]