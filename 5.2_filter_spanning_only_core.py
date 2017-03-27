import subprocess
import os
import argparse

parser = argparse.ArgumentParser('Filter sam/bam files to only keep spanning reads.')
parser.add_argument('--ref-genome', type = str, default = '/storage/resources/dbase/human/hs37d5/hs37d5.fa')
parser.add_argument('--out-pref', 	type = str, required = True)
parser.add_argument('--in-pref', 	type = str, required = True)
parser.add_argument('--locus-bed', 	type = str, required = True)
parser.add_argument('--read-len', 	type = int, required = True)


args = parser.parse_args()

ref_gen_dir = args.ref_genome
out_pref = args.out_pref
in_pref = args.in_pref
read_len = args.read_len
locus = args.locus_bed

# First extract locus start and end
with open(locus, 'r') as f:
	row = f.readline().split()
	chrom = row[0]
	locus_start = int(row[1])
	locus_end = int(row[2])

# Setting the filter edges for left and right
# Looking for reads such that r1 <= left_filter AND r2 >= right_filter
left_filter = locus_start - read_len
right_filter = locus_end

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
			if row[2] == chrom and int(row[3]) <= left_filter and int(row[7]) >= right_filter:
				out_sam_handle.write(record)

out_sam_handle.close()

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
