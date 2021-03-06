# This is the wrapper for the STR experiments.
# This script will take care of creating directories, and calls
# other scripts as necessary.
# Parsing and creating paths is dealt with in this script.

import errno    
import os
import sys
import subprocess

import argparse

# StackOverflow
def mkdir_p(path):
    try:
        os.makedirs(path)
        print 'Created path: ', path
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def find_motif_refAll(locus):
	with open(locus, 'r') as f:
		line = f.read().split()
		return line[4], (int(line[2]) - int(line[1]) + 1)

repo_dir = '/storage/nmmsv/str-expansions/'
base_dir = '/storage/nmmsv/expansion-experiments/'
run_dir = "/storage/nmmsv/analysis/GangSTR-analyses/simulation/"
ref = 'hg38'
if ref == 'hs37':
        ref_genome = '/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta'
        bed_dir = run_dir + "loci/hs37/" 
elif ref == 'hg38':
        ref_genome = '/storage/nmmsv/ref_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa'
        bed_dir = run_dir + "loci/hg38/" 
else:
        print 'Undefined reference genome.'


loci_file = '/storage/nmmsv/expansion-experiments/off_target_multi_locus_'+ref+'.txt'

loc_info = {}
with open (loci_file, 'r') as lfile:
        for lline in lfile:
                lrec = lline.strip().split('\t')
                loc = lrec[0]
                loc_info[loc] = {'chrom': lrec[1],\
                                 'start': int(lrec[2]),\
                                 'end': int(lrec[3]),\
                                 'motif_len': int(lrec[4]),\
                                 'motif': lrec[5],\
                                 'copy_list1': range(int(lrec[6]), int(lrec[7]), int(lrec[8])),\
                                 'copy_list2': range(int(lrec[9]), int(lrec[10]), int(lrec[11]))}

###################
coverage = 100
diploid = 'True'
number = 1
read_len = 150
#dist_mean  = range(read_len, 1100, read_len)
dist_mean = 500
dist_sdev  = 100

exp_name = 'offTarget_'+str(number) + \
           '_cov'+str(coverage)+\
           '_readLen_' + str(read_len)+\
           '_ref_' + ref

flank_len = 0
base_error = 0.005

mutat_rate = 0.005
indel_frac = 0.0001
indel_xtnd = 0.0001


num_threads = 5
bam_filter = True
heat_map_limit = 10


### Creating Directories ###
print '# Creating new directories..'
exp_dir = base_dir + exp_name + '/'
mkdir_p(base_dir)
mkdir_p(exp_dir)
############################


for lc in loc_info:
        locus = bed_dir+lc+'.bed'
        motif = loc_info[lc]['motif']
        str_len = loc_info[lc]['end'] - loc_info[lc]['start'] + 1
        ref_allele = int(str_len / len(motif))
        constant_allele = 0
        copy_list1 = loc_info[lc]['copy_list1']
        copy_list2 = loc_info[lc]['copy_list2']

        loc_dir = exp_dir + lc + '/'
        sim_gen_dir = loc_dir + '/simulated_genome/'
        sim_read_dir = loc_dir + '/simulated_read/'
        algn_read_dir = loc_dir + '/aligned_read/'
        temp_dir = loc_dir + '/temp/'
        mkdir_p(loc_dir)
        mkdir_p(sim_gen_dir)
        mkdir_p(sim_read_dir)
        mkdir_p(algn_read_dir)
        mkdir_p(temp_dir)


        ### STEP 0: Create Profile #####
        subprocess.call([	'python', 			repo_dir + '0_create_profile.py', \
                                '--exp-name',		exp_name, \
                                '--locus', 			locus, \
                                '--motif',			motif, \
                                '--flank-len',		str(flank_len), \
                                '--ref-gen-dir', 	ref_genome, \
                                '--repo-dir',		repo_dir, \
                                '--exp-dir',		exp_dir, \
                                '--read-len',		str(read_len), \
                                '--coverage',		str(coverage), \
                                '--read-ins-mean',	str(dist_mean), \
                                '--read-ins-stddev',str(dist_sdev), \
                                '--num-copy'] +		[str(nc) for nc in copy_list1] + \
                        ['--base-error',		str(base_error), \
                         '--num-threads',	str(num_threads), \
                         '--bam-filter', 	str(bam_filter), \
                         '--ref-allele-count',str(ref_allele), \
                         '--heat-map-limit', str(heat_map_limit), \
                         '--diploid',		diploid,\
                         '--constant-allele',str(constant_allele)])
        ################################


        ### STEP 1: Simulated Genome ###
        for nc1 in copy_list1:
                for nc2 in copy_list2:
                        out_path = sim_gen_dir + 'nc_' + str(nc1) + '_' + str(nc2) + '.fa'
                        subprocess.call(['python',repo_dir + '1.2_simulate_alt_genome_core.py', \
                                         '--ref-genome', ref_genome, \
                                         '--exp-name', 	exp_name, \
                                         '--out', 		out_path, \
                                         '--locus-bed',	locus, \
                                         '--motif',		motif, \
                                         '--flank-len', 	str(flank_len), \
                                         '--temp-dir', 	temp_dir, \
                                         '--num-copy',	str(nc1) ,\
                                         '--num-copy-2',	str(nc2) ,\
                                         '--diploid', 	diploid,\
                                         '--exp-dir',	exp_dir])
                                
        #############################

        ### STEP 2: wgsim ###########
        for nc1 in copy_list1:
                for nc2 in copy_list2:
                        in_path = sim_gen_dir + 'nc_' + str(nc1) + '_' + str(nc2) + '.fa'
                        out_pref= sim_read_dir + 'nc' + str(nc1) + '_' + str(nc2)
                        subprocess.call(['python',repo_dir + '2.2_read_simulated_data_core.py', \
                                         '--exp-name', 	exp_name, \
                                         '--fasta-in',	in_path, \
                                         '--out-pref', 	out_pref, \
                                         '--coverage',	str(coverage), \
                                         '--num-copy',	str(nc), \
                                         '--dist-mean', 	str(dist_mean), \
                                         '--dist-sdev', 	str(dist_sdev), \
                                         '--motif',		motif, \
                                         '--base-error',	str(base_error), \
                                         '--flank-len', 	str(flank_len), \
                                         '--read-len',	str(read_len), \
                                         '--mutat-rate',	str(mutat_rate), \
                                         '--indel-frac', str(indel_frac), \
                                         '--indel-xtnd', str(indel_xtnd)])

        #############################

        ### STEP 3: Alignment #######
        with open(algn_read_dir + 'bamlist.txt', 'w') as bamlist:
                for nc1 in copy_list1:
                        for nc2 in copy_list2:
                                read_grp_header = 	'\'@RG\\tID:' + exp_name + \
                                                        '\\tSM:' + str(nc1) + '_' + str(nc2) + \
                                                        '\\tLB:' + str(coverage)+ \
                                                        '\\tPL:' + str(base_error) + '\''
                                in_pref = sim_read_dir + 'nc' + str(nc1) + '_' + str(nc2)
                                out_pref = algn_read_dir + 'nc' + str(nc1) + '_' + str(nc2)
                                subprocess.call(['python',repo_dir + '3.2_align_read_core.py', \
                                                 '--ref-genome', ref_genome, \
                                                 '--out-pref', 	out_pref, \
                                                 '--in-pref', 	in_pref, \
                                                 '--read-grp',	read_grp_header, \
                                                 '--num-threads',str(num_threads), \
                                                 '--cpp-data',	str(1)])
                                bamlist.write(out_pref+'.sorted.bam\n')

        ##############################
