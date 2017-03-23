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
parser.add_argument('--read-len', type = int, default = 100)
parser.add_argument('--num-reads', type = int, nargs = '+', default = [1000])
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
expName = args.exp_name
repoDir = args.repo_dir
mkdir_p(repoDir + 'experiments/')
expCaseDir = repoDir + 'experiments/' + expName
mkdir_p(expCaseDir)

with open(expCaseDir + '/profile.txt', 'w') as f:
	pickle.dump(vars(args), f)

# if len(sys.argv) < 2:
# 	print '### Usage python 0_create_profile.py dir expName locus motif flankLength refGenomeDir nCopyCount nCopy1 nCopy2 ...'
# 	sys.exit()
# else:
# 	repoDir = sys.argv[1]
# 	expName = sys.argv[2]
# 	expCaseDir = repoDir + 'experiments/' + expName
# 	mkdir_p(repoDir + 'experiments/' + expCaseDir)
# 	expProfile = {'motif': sys.argv[3],\
# 				'locus': 'loci/' + sys.argv[4] + '.bed',\
# 				'flankLength': int(sys.argv[5]),\
# 				'refGenomeDir': sys.argv[6]}
# 	motif = 
# 	locus = 'loci/' + sys.argv[4] + '.bed'
# 	flankLength = int(sys.argv[5])
# 	refGenomeDir = sys.argv[6]
# 	nCopyCount = int(sys.argv[7])
# 	nCopyList = []
# 	for i in range(8, 8 + nCopyCount):
# 		nCopyList.append(int(sys.argv[i]))

# 	with (expCaseDir + '/exp_profile.txt', 'w') as f:
# 		f.write('#### Experiment Profile ####')
# 		f.write('>>locus:\t', locus)
# 		f.write('>>motif:\t', motif)
# 		f.write('>>flankLength:\t', flankLength)
