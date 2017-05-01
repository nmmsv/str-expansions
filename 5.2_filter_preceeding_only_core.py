import subprocess
import os, sys
import argparse
sys.path.append('/storage/nmmsv/str-expansions/functions/')
from realignment import smith_waterman
from load_info import load_profile, extract_locus_info

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
# First extract locus start and end
# with open(locus, 'r') as f:
# 	row = f.readline().split()
# 	chrom = row[0]
# 	locus_start = int(row[1])
# 	locus_end = int(row[2])

chrom, locus_start, locus_end = extract_locus_info(locus)

score_dict = {	'match': 	3, \
				'mismatch': -1, \
				'gap': 		-3}
verbose = False


# Setting the filter edges for left and right
# Looking for reads such that r1 <= left_filter AND r2 >= right_filter
left_filter = locus_start - read_len
right_filter = locus_end

realign = realignment_string (temp_fa_dir, read_len, motif)

in_sam = in_pref + '.sam'
out_sam = out_pref + '.sam'
out_sam_handle = open(out_sam, 'w')
print 'Filtering ' + in_pref + '.sam'

list_rows = []
with open(in_sam, 'r') as in_sam_handle:
	for record in in_sam_handle:
		if record[0] == '@':
			out_sam_handle.write(record)
		else:
			row = record.split()
			if (row[2] == chrom and row[6] == '=') or (row[6] == chrom):		# checking correct chrom for mate 
				if int(row[7]) <= locus_start:									# checking if mate is before STR region
					if row[2] != chrom:
						out_sam_handle.write(record)
						print row[9]
						pos = smith_waterman(realign, row[9], score_dict, verbose)
						if pos > read_len - 2 and pos < read_len + 4:
							print '## IRR'
						elif pos <= read_len - 2:
							print '## Pre-Flank', row[0]
						else:
							print '## Post-Flank (mate before STR)', row[0]
					elif int(row[3]) >= locus_start - read_len and int(row[3]) <= locus_end:
						# list_rows.append(row[0])

						out_sam_handle.write(record)
						print row[9]
						pos = smith_waterman(realign, row[9], score_dict, verbose)
						if pos > read_len - 2 and pos < read_len + 4:
							print '## IRR'
						elif pos <= read_len - 2:
							print '## Pre-Flank', row[0]
						else:
							print '## Post-Flank (mate before STR)', row[0]


out_sam_handle.close()
print locus_start - read_len
print locus_end
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
