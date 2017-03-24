import os
import argparse

parser = argparse.ArgumentParser('Align simulated reads using bwa mem.')
parser.add_argument('--ref-genome', type = str, default = '/storage/resources/dbase/human/hs37d5/hs37d5.fa')
parser.add_argument('--exp-name', 	type = str, required = True)
parser.add_argument('--out-pref', 	type = str, required = True)
parser.add_argument('--in-pref', 	type = str, required = True)
parser.add_argument('--read-grp', 	type = str, required = True)
parser.add_argument('--num-threads',type = int)

args = parser.parse_args()

exp_name = args.exp_name
ref_gen_dir = args.ref_genome
out_pref = args.out_pref
in1_file = args.in_pref + '.read1.fq'
in2_file = args.in_pref + '.read2.fq'
num_thrd = args.num_threads
read_grp = args.read_grp

#subprocess.call(['bwa', 'mem',\
#		refGenomeDir, inFq1, inFq2,\
#		'-R', '@RG\\tID:HTT\\tSM:16\\tLB:lb\\tPL:pl',\
#		'>', outFile + '.sam'])
print '# Executing bwa mem on ' + args.in_pref
os.system('bwa mem ' + \
		ref_gen_dir + ' ' +  \
		in1_file + ' ' + \
		in2_file + ' ' + \
		'-t ' + str(num_thrd) + ' ' + \
		'-R ' + str(read_grp) + ' ' + \
		'> ' + out_pref + '.sam')
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
