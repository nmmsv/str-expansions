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
parser.add_argument('--loci-beds', 	type = str, nargs = '+', required = True)
parser.add_argument('--motifs', 	type = str, nargs = '+', required = True)
parser.add_argument('--ref-alleles', 	type = int, nargs = '+', required = True)
parser.add_argument('--flank-len', 	type = int, required = True)
parser.add_argument('--temp-dir', 	type = str, required = True)
parser.add_argument('--num-copy', 	type = int, required = True)
parser.add_argument('--num-copy-2',	type = int, required = False)
parser.add_argument('--diploid', 	type = str, required = True)
parser.add_argument('--exp-dir', 	type = str, required = True)
args = parser.parse_args()

loci_beds = args.loci_beds
ref_gen_dir = args.ref_genome
ref_alleles = args.ref_alleles
motifs = args.motifs
flank_len = args.flank_len
num_copy = args.num_copy
out_file = args.out
exp_name = args.exp_name
temp_dir = args.temp_dir
diploid = args.diploid
exp_dir = args.exp_dir


arg_dict = load_profile(exp_dir)
if args.num_copy_2 is not None:
	constant_allele = args.num_copy_2
else:
	constant_allele = arg_dict['constant_allele']
#ref_allele = arg_dict['ref_allele_count']
# STR_locus contains the information for the STR locus
# Format: Chr	Start	End	MotifLength	nCopies
loci_names = [loc.split('/')[-1].split('.')[0] for loc in loci_beds]

chroms = []
start_locs = []
end_locs = []
begin_flanks = []
end_flanks = []
for i in range(len(loci_beds)):
        locus_bed = loci_beds[i]
        with open (locus_bed, 'r') as f:
                for row in csv.reader(f, dialect='excel-tab'):
                        print 'STR locus: ', row
                        chroms.append(row[0])
                        start_locs.append(int(row[1]))
                        end_locs.append(int(row[2]))
                        begin_flanks.append(int(row[1]) - flank_len)
                        end_flanks.append(int(row[2]) + flank_len)
                        break



print '# Creating temp bed files..'
refHaps = []
prefix_haplos = []
suffix_haplos = []
for i in range(len(begin_flanks)):
        ref_bed = open(temp_dir + '/tempRefBed'+str(i)+'.bed', 'w')
        ref_bed.write(chroms[i] + '\t' + str(begin_flanks[i]) + '\t' + str(end_flanks[i]) + '\n')
        ref_bed.close()

        prefix_bed = open(temp_dir + '/tempPrefixBed'+str(i)+'.bed', 'w')
        prefix_bed.write(chroms[i] + '\t' + str(begin_flanks[i]) + '\t' + str(start_locs[i] - 1) + '\n')
        prefix_bed.close()

        suffix_bed = open(temp_dir + '/tempSuffixBed'+str(i)+'.bed', 'w')
        suffix_bed.write(chroms[i] + '\t' + str(end_locs[i]) + '\t' + str(end_flanks[i]) + '\n')
        suffix_bed.close()


        ref_allele_fa = temp_dir + exp_name + '_ref'+str(i)+'.fa'
        prefix_fa = temp_dir + exp_name + '_prefix'+str(i)+'.fa'
        suffix_fa = temp_dir + exp_name + '_suffix'+str(i)+'.fa'
        # Create blank fasta file if it doesn't exist:
        open(ref_allele_fa, 'a').close()
        open(prefix_fa, 'a').close()
        open(suffix_fa, 'a').close()
        # Extracting reference haplotype
        subprocess.call(['bedtools', 'getfasta', \
                         '-fi', ref_gen_dir, \
                         '-bed', temp_dir + '/tempRefBed'+str(i)+'.bed', \
                         '-fo', ref_allele_fa])
        with open (ref_allele_fa, 'r') as refFasta:
                refFasta.readline().strip()
                refHaps.append(refFasta.readline().strip())

        print '# Calling bedtools getfasta..'
        # Generating prefix and suffix for altered (STR nCopy changed) haplotype
        subprocess.call(['bedtools', 'getfasta', \
                         '-fi', ref_gen_dir, \
                         '-bed', temp_dir + 'tempPrefixBed'+str(i)+'.bed', \
                         '-fo', prefix_fa])
        with open (prefix_fa, 'r') as preFile:
                preFile.readline().strip()
                prefix_haplos.append(preFile.readline().strip())

        subprocess.call(['bedtools', 'getfasta', \
                         '-fi', ref_gen_dir, \
                         '-bed', temp_dir + 'tempSuffixBed'+str(i)+'.bed', \
                         '-fo', suffix_fa])
        with open (suffix_fa, 'r') as sufFile:
                sufFile.readline().strip()
                suffix_haplos.append(sufFile.readline().strip())

	
with open(out_file, 'w') as f:
        for i in range(len(ref_alleles)):
                ref_allele = ref_alleles[i]
                if diploid == 'True' and constant_allele == ref_allele:
                        f.write('>' + exp_name + '_' + str(num_copy) + '_altAllele_'+loci_names[i]+'\n')
                        f.write(prefix_haplos[i] + motifs[i] * num_copy + suffix_haplos[i] + '\n')
                        f.write('>' + exp_name + '_' + str(num_copy)+'_refAllele_'+loci_names[i]+'\n')
                        f.write(refHaps[i] + '\n')
                elif diploid == 'True' and constant_allele != ref_allele:
                        f.write('>' + exp_name + '_' + str(num_copy) + '_altAllele_'+loci_names[i]+'\n')
                        f.write(prefix_haplos[i] + motifs[i] * num_copy + suffix_haplos[i] + '\n')
                        f.write('>' + exp_name + '_' + str(constant_allele)+'_constAllele_'+loci_names[i]+'\n')
                        f.write(prefix_haplos[i] + motifs[i] * constant_allele + suffix_haplos[i] + '\n')
                else:
                        f.write('>' + exp_name + '_' + str(num_copy) + '_haplo_'+loci_names[i]+'\n')
                        f.write(prefix_haplos[i] + motifs[i] * num_copy + suffix_haplos[i] + '\n')
print '# Created fasta file: ' + out_file


