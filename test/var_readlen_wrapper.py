# This is the wrapper for the STR experiments.
# This script will take care of creating directories, and calls
# other scripts as necessary.
# Parsing and creating paths is dealt with in this script.

import errno    
import os
import sys
import subprocess

import argparse

def find_motif_refAll(locus):
	with open(locus, 'r') as f:
		line = f.read().split()
		return line[4], (int(line[2]) - int(line[1]) + 1)

# parser = argparse.ArgumentParser('Filter sam/bam files to only keep spanning reads.')
# # parser.add_argument('--align', type = str, required = True)

# args = parser.parse_args()

align_flag = 'True'
# if align_flag == 'False':
# 	print ""
# 	print "####### NO ALIGNMENT!!!!"
# 	print ""
# 	print "########################"
# 	print "########################"
# Parameters ######

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


###################
locus_name = 'C9ORF72'
coverage = 40
diploid = 'True'
number = 4
read_len = range(50,305,50)
#dist_mean  = range(read_len, 1100, read_len)
dist_mean = 500
dist_sdev  = 100
#copy_list1 = range(10,100,12) + range(110, 500, 20)
#copy_list2 =  [20, 30, 40]
copy_list1 = range(5,1010,100)
copy_list2 = range(5, 31, 5)


exp_name = 'var_readlen_'+str(number)+'_' + locus_name + \
			'_cov'+str(coverage)+\
			'_distmean_' + str(dist_mean)


###################

locus = bed_dir+locus_name+'.bed'
motif, str_len = find_motif_refAll(locus)

ref_allele = int(str_len / len(motif))
constant_allele = 0
flank_len = 10000
base_error = 0.005

mutat_rate = 0.0
indel_frac = 0.0
indel_xtnd = 0.0


num_threads = 5
bam_filter = True
heat_map_limit = 10
###################


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

### Creating Directories ###
print '# Creating new directories..'
exp_dir = base_dir + exp_name + '/'
sim_gen_dir = exp_dir + '/simulated_genome/'
sim_read_dir = exp_dir + '/simulated_read/'
algn_read_dir = exp_dir + '/aligned_read/'
# estm_dir = exp_dir + '/estimation/'
# plot_dir = exp_dir + '/plots/'
temp_dir = exp_dir + '/temp/'
mkdir_p(base_dir)
mkdir_p(exp_dir)
mkdir_p(sim_gen_dir)
mkdir_p(sim_read_dir)
mkdir_p(algn_read_dir)
# mkdir_p(estm_dir)
# mkdir_p(plot_dir)
mkdir_p(temp_dir)
############################


### STEP 0: Create Profile #####
if align_flag == 'True':
	subprocess.call([	'python', 			repo_dir + '0_create_profile.py', \
						'--exp-name',		exp_name, \
						'--locus', 			locus, \
						'--motif',			motif, \
						'--flank-len',		str(flank_len), \
						'--ref-gen-dir', 	ref_genome, \
						'--repo-dir',		repo_dir, \
						'--exp-dir',		exp_dir, \
						'--read-len',		str(max(read_len)), \
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
if align_flag == 'True':
	for nc1 in copy_list1:
		for nc2 in copy_list2:
			out_path = sim_gen_dir + 'nc_' + str(nc1) + '_' + str(nc2) + '.fa'
			subprocess.call(['python', 		repo_dir + '1.2_simulate_alt_genome_core.py', \
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
if align_flag == 'True':
        for rl_index in range(len(read_len)):
                for nc1 in copy_list1:
                        for nc2 in copy_list2:
                                in_path = sim_gen_dir + 'nc_' + str(nc1) + '_' + str(nc2) + '.fa'
                                out_pref= sim_read_dir + 'is' + str(read_len[rl_index]) + '_nc' + str(nc1) + '_' + str(nc2)
                                subprocess.call(['python', 		repo_dir + '2.2_read_simulated_data_core.py', \
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
                                                 '--read-len',	str(read_len[rl_index]), \
                                                 '--mutat-rate',	str(mutat_rate), \
                                                 '--indel-frac', str(indel_frac), \
                                                 '--indel-xtnd', str(indel_xtnd)])

#############################

### STEP 3: Alignment #######
with open(algn_read_dir + 'readlen_list.txt', 'w') as islist:
        for rl_index in range(len(read_len)):
                islist.write(str(read_len[rl_index]) + '\n')
                with open(algn_read_dir + str(read_len[rl_index]) + '_bamlist.txt', 'w') as bamlist:
                        for nc1 in copy_list1:
                                for nc2 in copy_list2:
                                        read_grp_header = 	'\'@RG\\tID:' + exp_name + \
                                                                '\\tSM:' + str(nc1) + '_' + str(nc2) + \
                                                                '\\tLB:' + str(coverage)+ \
                                                                '\\tPL:' + str(base_error) + '\''
                                        in_pref = sim_read_dir + 'is' + str(read_len[rl_index]) + '_nc' + str(nc1) + '_' + str(nc2)
                                        out_pref = algn_read_dir + 'is' + str(read_len[rl_index]) + '_nc' + str(nc1) + '_' + str(nc2)
                                        subprocess.call(['python', 		repo_dir + '3.2_align_read_core.py', \
                                                         '--ref-genome', ref_genome, \
                                                         '--out-pref', 	out_pref, \
                                                         '--in-pref', 	in_pref, \
                                                         '--read-grp',	read_grp_header, \
                                                         '--num-threads',str(num_threads), \
                                                         '--cpp-data',	str(1)])
                                        bamlist.write(out_pref+'.sorted.bam\n')

##############################


# ### STEP 5: Filter: Find Spanning Read Pairs #######
# for nc in copy_list:
# 	in_pref = algn_read_dir + 'nc_' + str(nc)
# 	out_pref = algn_read_dir + 'nc_' + str(nc) + '_srp'
# 	subprocess.call(['python', 		repo_dir + '5.2_filter_spanning_only_core.py', \
# 					'--out-pref', 	out_pref, \
# 					'--in-pref', 	in_pref, \
# 					'--exp-dir',	exp_dir ])

# ##################################################################

# ### STEP 5: Filter: Find Enclosing Reads #######
# for nc in copy_list:
# 	in_pref = algn_read_dir + 'nc_' + str(nc)
# 	out_pref = algn_read_dir + 'nc_' + str(nc) + '_er'
# 	subprocess.call(['python', 		repo_dir + '5.2_filter_enclosing_only_core.py', \
# 					'--out-pref', 	out_pref, \
# 					'--in-pref', 	in_pref, \
# 					'--exp-dir',	exp_dir ])

# #################################################################

# ### STEP 5: Filter: Find Fully Repetitive Reads #######
# for nc in copy_list:
# 	in_pref = algn_read_dir + 'nc_' + str(nc)
# 	out_pref = algn_read_dir + 'nc_' + str(nc) + '_frr'
# 	subprocess.call(['python', 		repo_dir + '5.2_filter_FRR_only_core.py', \
# 					'--out-pref', 	out_pref, \
# 					'--in-pref', 	in_pref, \
# 					'--exp-dir',	exp_dir ])

##################################################################

# ### STEP 9: Estimate: Estimate STR length, compute error, and plot #####
# for nc in copy_list:
# 	in_path = algn_read_dir + 'nc_' + str(nc) + '_flt.bam'
# 	plot_path = plot_dir + 'nc_' + str(nc) + '.pdf'
# 	estm_path = estm_dir + 'nc_' + str(nc) + '.txt'
# 	subprocess.call(['python', 		repo_dir + '9_estimate_str_and_plot_core.py', \
# 					'--in-path', 	in_path, \
# 					'--estm-path',	estm_path,\
# 					'--plot-path',	plot_path,\
# 					'--ref-allele',	str(ref_allele),\
# 					'--alt-allele',	str(nc),\
# 					'--motif', 		motif,\
# 					'--dist-mean',	str(dist_mean),\
# 					'--temp-dir',	temp_dir])

# #########################################################################
