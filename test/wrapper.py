# This is the wrapper for the STR experiments.
# This script will take care of creating directories, and calls
# other scripts as necessary.
# Parsing and creating paths is dealt with in this script.

import errno    
import os
import sys
import subprocess


# Parameters ######

repo_dir = '/storage/nmmsv/str-expansions/'
base_dir = '/storage/nmmsv/expansion-experiments/'
ref_genome = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
locus = repo_dir + '/loci/ATXN7.bed'
motif = 'GCA'


exp_name = 'ATXN7_12_core'
nCopyList = [3, 10, 18]
flank_len = 1000
base_error = 0.0
dist_mean  = 1000
dist_sdev  = 100
coverage = 100
read_len = 100
mutat_rate = 0.0
indel_frac = 0.0
indel_xtnd = 0.0

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
temp_dir = exp_dir + '/temp/'
mkdir_p(base_dir)
mkdir_p(exp_dir)
mkdir_p(sim_gen_dir)
mkdir_p(sim_read_dir)
mkdir_p(algn_read_dir)
mkdir_p(temp_dir)
############################

### STEP 1: Simulated Genome ###

for nc in nCopyList:
	out_path = sim_gen_dir + 'nc_' + str(nc) + '.fa'
	subprocess.call(['python', 		repo_dir + '1.2_simulate_alt_genome_core.py', \
					'--ref-genome', '/storage/resources/dbase/human/hs37d5/hs37d5.fa', \
					'--exp-name', 	exp_name, \
					'--out', 		out_path, \
					'--locus-bed',	locus, \
					'--motif',		motif, \
					'--flank-len', 	str(flank_len), \
					'--temp-dir', 	temp_dir, \
					'--num-copy',	str(nc)])

#############################

### STEP 2: wgsim ###########
for nc in nCopyList:
	in_path = sim_gen_dir + 'nc_' + str(nc) + '.fa'
	out_pref= sim_read_dir + 'nc_' + str(nc)
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
					'--read-len',	str(read_len), \
					'--mutat-rate',	str(mutat_rate), \
					'--indel-frac', str(indel_frac), \
					'--indel-xtnd', str(indel_xtnd)])

#############################

### STEP 3: Alignment #######

for nc in nCopyList:
	read_grp_header = 	'\'@RG\\tID:' + exp_name + \
					'\\tSM:' + str(nc) + \
					'\\tLB:' + str(coverage)+ \
					'\\tPL:' + str(base_error) + '\''
	in_pref = sim_read_dir + 'nc_' + str(nc)
	out_pref = algn_read_dir + 'nc_' + str(nc)
	subprocess.call(['python', 		repo_dir + '3.2_align_read_core.py', \
					'--ref-genome', '/storage/resources/dbase/human/hs37d5/hs37d5.fa', \
					'--exp-name', 	exp_name, \
					'--out-pref', 	out_pref, \
					'--in-pref', 	in_pref, \
					'--read-grp',	read_grp_header, \
					'--num-threads',str(5)])

##############################


