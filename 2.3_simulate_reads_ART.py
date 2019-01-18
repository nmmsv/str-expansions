import subprocess
import sys
import numpy as np
import argparse


parser = argparse.ArgumentParser('Generate simulated reads using wgsim.')
parser.add_argument('--exp-name', 	type = str, 	required = True)
parser.add_argument('--fasta-in', 	type = str, 	required = True)
parser.add_argument('--out-pref', 	type = str, 	required = True)
parser.add_argument('--motif', 		type = str, 	required = True)
parser.add_argument('--base-error', type = float, 	required = True)
parser.add_argument('--dist-mean', 	type = int, 	required = True)
parser.add_argument('--dist-sdev', 	type = int, 	required = True)
parser.add_argument('--flank-len', 	type = int, 	required = True)
parser.add_argument('--coverage', 	type = int, 	required = True)
parser.add_argument('--read-len', 	type = int, 	required = True)
parser.add_argument('--num-copy', 	type = int, 	required = True)
parser.add_argument('--mutat-rate', type = float, 	required = True)
parser.add_argument('--indel-frac', type = float, 	required = True)
parser.add_argument('--indel-xtnd', type = float, 	required = True)
args = parser.parse_args()

exp_name 	= args.exp_name
in_file		= args.fasta_in
out_pref 	= args.out_pref + '.'
motif 		= args.motif
base_error 	= args.base_error
dist_mean	= args.dist_mean
dist_sdev	= args.dist_sdev
flank_len 	= args.flank_len
coverage 	= args.coverage
read_len 	= args.read_len
num_copy 	= args.num_copy
mutat_rate	= args.mutat_rate
indel_frac 	= args.indel_frac
indel_xtnd	= args.indel_xtnd

num_reads = coverage * (2 * flank_len + num_copy * len(motif)) / 2 / read_len
seed = int(np.random.rand(1)[0] * 100)
#seed = 60

# HiSeqX PCR free (150bp)
subprocess.call(['/storage/nmmsv/ART/art_bin_MountRainier/'+\
                 'art_illumina',\
                 '-ss', 'HSXn',\
                 '-i', in_file,\
                 '-l', str(read_len),\
                 '-f', str(coverage),\
                 '-o', out_pref,\
                 '-m', str(dist_mean),\
                 '-s', str(dist_sdev),\
                 '-na', '-p'])


