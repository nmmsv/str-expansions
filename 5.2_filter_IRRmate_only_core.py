import subprocess
import os, sys
import argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import expansion_aware_realign, classify_realigned_read
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


chrom, locus_start, locus_end = extract_locus_info(locus)
pre, post = extract_pre_post_flank(exp_dir, read_len)

score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False
margin = 2

# Setting the filter edges for left and right
# Looking for reads such that r1 <= left_filter AND r2 >= right_filter
left_filter = locus_start - read_len
right_filter = locus_end

realign = realignment_string (temp_fa_dir, read_len, motif)

in_sam = in_pref + '.sam'
print 'Filtering ' + in_pref + '.sam'

list_reads = []
with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] != '@':
			row = record.split()
			if len(row[9]) > 0.9 * read_len:		# Discarding soft clipped reads
				if (row[2] == chrom and row[6] == '=') or (row[6] == chrom):		# checking correct chrom for mate 
					if int(row[7]) <= locus_start - read_len:						# checking if mate is before STR region
						if row[2] != chrom:
							# Performing realignment to check if IRR
							sample = row[9]
							nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
							read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
							if read_class == 'IRR':
								list_reads.append(row[0])			# Saving read_ID of potential reads that going to be extracted
								print 'Found IRR!'
							else:
								pass
						elif int(row[3]) >= locus_start and int(row[3]) <= locus_end - read_len:
							# Performing realignment to check if IRR
							sample = row[9]
							nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
							read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
							if read_class == 'IRR':
								list_reads.append(row[0])			# Saving read_ID of potential reads that going to be extracted
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
			if row[0] in list_reads:
				# perform realignment to check this read is the non-spanning read (not the IRR mate)
				sample = row[9]
				nCopy, pos, score = expansion_aware_realign(sample, pre, post, motif, score_dict, verbose)
				read_class = classify_realigned_read(sample, motif, pos, nCopy, score, score_dict, read_len, margin, verbose)
				if read_class == 'NoSpan':
					out_sam_handle.write(record)
					print 'Record Added!'



out_sam_handle.close()

# print list_rows
# Sort and index filtered bam
os.system('samtools view -bT ' + \
		ref_gen_dir + ' ' + \
		out_pref + '.sam' + ' ' + 
		'> ' + out_pref + '.bam')
os.system('samtools sort -o ' + \
		out_pref + '.sorted.bam '+ \
		out_pref + '.bam')
os.system('samtools index ' + \
		out_pref + '.sorted.bam ' + \
		out_pref + '.sorted.bai')
os.system('samtools index ' + out_pref + '.sorted.bam')
