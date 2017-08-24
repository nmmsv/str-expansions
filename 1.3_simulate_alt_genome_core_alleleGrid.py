import csv
import subprocess
import pickle
import sys
import argparse

sys.path.append('/storage/nmmsv/str-expansions/functions/')
from load_info import load_profile, extract_locus_info

parser = argparse.ArgumentParser('Simulate alterred genome using bedtools getfasta.')
parser.add_argument('--ref-genome', type = str, default = '/storage/resources/dbase/human/hs37d5/hs37d5.fa')
parser.add_argument('--exp-name', 	type = str, required = True)
parser.add_argument('--out', 		type = str, required = True)
parser.add_argument('--locus-bed', 	type = str, required = True)
parser.add_argument('--motif', 		type = str, required = True)
parser.add_argument('--flank-len', 	type = int, required = True)
parser.add_argument('--temp-dir', 	type = str, required = True)
parser.add_argument('--num-copy', 	type = int, required = True)
parser.add_argument('--diploid', 	type = str, required = True)
parser.add_argument('--exp-dir', 	type = str, required = True)
args = parser.parse_args()

locus_bed = args.locus_bed
ref_gen_dir = args.ref_genome
motif = args.motif
flank_len = args.flank_len
num_copy = args.num_copy
out_file = args.out
exp_name = args.exp_name
temp_dir = args.temp_dir
diploid = args.diploid
exp_dir = args.exp_dir

arg_dict = load_profile(exp_dir)
constant_allele = arg_dict['constant_allele']
ref_allele = arg_dict['ref_allele_count']
# STR_locus contains the information for the STR locus
# Format: Chr	Start	End	MotifLength	nCopies
with open (locus_bed, 'r') as f:
	for row in csv.reader(f, dialect='excel-tab'):
		print 'STR locus: ', row
		chrom = row[0]
		start_loc = int(row[1])
		end_loc = int(row[2])


begin_flank = start_loc - flank_len
end_flank = end_loc + flank_len

print '# Creating temp bed files..'
ref_bed = open(temp_dir + '/tempRefBed.bed', 'w')
ref_bed.write(chrom + '\t' + str(begin_flank) + '\t' + str(end_flank) + '\n')
ref_bed.close()

prefix_bed = open(temp_dir + 'tempPrefixBed.bed', 'w')
prefix_bed.write(chrom + '\t' + str(begin_flank) + '\t' + str(start_loc - 1) + '\n')
prefix_bed.close()

suffix_bed = open(temp_dir + 'tempSuffixBed.bed', 'w')
suffix_bed.write(chrom + '\t' + str(end_loc) + '\t' + str(end_flank) + '\n')
suffix_bed.close()

ref_allele_fa = temp_dir + exp_name + '_ref.fa'
prefix_fa = temp_dir + exp_name + '_prefix.fa'
suffix_fa = temp_dir + exp_name + '_suffix.fa'
# Create blank fasta file if it doesn't exist:
open(ref_allele_fa, 'a').close()
open(prefix_fa, 'a').close()
open(suffix_fa, 'a').close()
# Extracting reference haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', ref_gen_dir, \
	'-bed', temp_dir + '/tempRefBed.bed', \
	'-fo', ref_allele_fa])
with open (ref_allele_fa, 'r') as refFasta:
	refFasta.readline().strip()
	refHap = refFasta.readline().strip()

print '# Calling bedtools getfasta..'
# Generating prefix and suffix for altered (STR nCopy changed) haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', ref_gen_dir, \
	'-bed', temp_dir + 'tempPrefixBed.bed', \
	'-fo', prefix_fa])
with open (prefix_fa, 'r') as preFile:
	preFile.readline().strip()
	prefix_haplo = preFile.readline().strip()

subprocess.call(['bedtools', 'getfasta', \
	'-fi', ref_gen_dir, \
	'-bed', temp_dir + 'tempSuffixBed.bed', \
	'-fo', suffix_fa])
with open (suffix_fa, 'r') as sufFile:
	sufFile.readline().strip()
	suffix_haplo = sufFile.readline().strip()

	
with open(out_file, 'w') as f:
	if diploid == 'True' and constant_allele == ref_allele:
		f.write('>' + exp_name + '_' + str(num_copy) + '_altAllele\n')
		f.write(prefix_haplo + motif * num_copy + suffix_haplo + '\n')
		f.write('>' + exp_name + '_' + str(num_copy)+'_refAllele\n')
		f.write(refHap + '\n')
	elif diploid == 'True' and constant_allele != ref_allele:
		f.write('>' + exp_name + '_' + str(num_copy) + '_altAllele\n')
		f.write(prefix_haplo + motif * num_copy + suffix_haplo + '\n')
		f.write('>' + exp_name + '_' + str(constant_allele)+'_constAllele\n')
		f.write(prefix_haplo + motif * constant_allele + suffix_haplo + '\n')
	else:
		f.write('>' + exp_name + '_' + str(num_copy) + '_haplo\n')
		f.write(prefix_haplo + motif * num_copy + suffix_haplo + '\n')
print '# Created fasta file: ' + out_file


