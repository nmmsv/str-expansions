import os
import argparse

parser = argparse.ArgumentParser('Align simulated reads using bwa mem.')
parser.add_argument('--ref-genome', type = str, default = '/storage/resources/dbase/human/hs37d5/hs37d5.fa')
parser.add_argument('--out-pref', 	type = str, required = True)
parser.add_argument('--in-pref', 	type = str, required = True)
parser.add_argument('--read-grp', 	type = str, required = True)
parser.add_argument('--num-threads',type = int)
parser.add_argument('--cpp-data', 	type = int, default = 0)

args = parser.parse_args()

ref_gen_dir = args.ref_genome
out_pref = args.out_pref
in1_file = args.in_pref + '.1.fq'
in2_file = args.in_pref + '.2.fq'
num_thrd = args.num_threads
read_grp = args.read_grp
cpp_data = int(args.cpp_data)

print '# Executing bwa mem on ' + args.in_pref

print
print ref_gen_dir

os.system('bwa mem ' + \
		ref_gen_dir + ' ' +  \
		in1_file + ' ' + \
		in2_file + ' ' + \
          '-M ' + \
		'-t ' + str(num_thrd) + ' ' + \
		'-R ' + str(read_grp) + ' ' + \
		'> ' + out_pref + '.sam')
if cpp_data == 1:
	os.system('samtools view -bT ' + \
			ref_gen_dir + ' ' + \
			out_pref + '.sam' + ' ' + 
			'> ' + out_pref + '.bam')
	os.system('samtools sort -o ' + \
			out_pref + '.sorted.bam '+ \
			out_pref + '.bam')
	os.system('samtools index ' + \
			out_pref + '.sorted.bam ')
	os.system('rm ' + out_pref + '.sam')
	os.system('rm ' + out_pref + '.bam')
# os.system('samtools index ' + out_pref + '.sorted.bam')
