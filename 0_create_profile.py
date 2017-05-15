import sys
import errno    
import os
import sys
import pickle
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

parser = argparse.ArgumentParser('Creating profile for str-expansions experiments')
parser.add_argument('--exp-name', type = str, required = True)
parser.add_argument('--locus', type = str, required = True)
parser.add_argument('--motif', type = str, required = True)
parser.add_argument('--flank-len', type = int, default = 1000)
parser.add_argument('--ref-gen-dir', type = str, \
					default = '/storage/resources/dbase/human/hs37d5/hs37d5.fa')
parser.add_argument('--repo-dir', type = str, \
					default = '/storage/nmmsv/str-expansions/')
parser.add_argument('--exp-dir', type = str, \
                    default = '/strorage/expansion-experiments/')
parser.add_argument('--read-len', type = int, default = 100)
parser.add_argument('--num-reads', type = int, nargs = '+', default = [1000]) # deprecated
# Modified the line bellow to fit heatmap needs
parser.add_argument('--coverage', type = int, nargs = '+', required = True)
parser.add_argument('--read-ins-mean', type = int, default = 500)
parser.add_argument('--read-ins-stddev', type = int, default = 50)
parser.add_argument('--num-copy', type = int, nargs = '+', required = True)
parser.add_argument('--base-error', type = float, default = 0.0)
parser.add_argument('--num-threads', type = int, default = 4)
parser.add_argument('--bam-filter', type = bool, default = True)
parser.add_argument('--ref-allele-count', type = int, required = True)
parser.add_argument('--heat-map-limit', type = int, default = 10)
args = parser.parse_args()

print vars(args)
exp_dir = args.exp_dir

with open(exp_dir + '/profile.txt', 'w') as f:
	pickle.dump(vars(args), f)

