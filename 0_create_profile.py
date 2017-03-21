import sys
import errno    
import os
import sys
import pickle
import argeparse


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

parser = argeparse.ArgumentParser('Creating profile for ')


if len(sys.argv) < 2:
	print '### Usage python 0_create_profile.py dir expName locus motif flankLength refGenomeDir nCopyCount nCopy1 nCopy2 ...'
	sys.exit()
else:
	repoDir = sys.argv[1]
	expName = sys.argv[2]
	expCaseDir = repoDir + 'experiments/' + expName
	mkdir_p(repoDir + 'experiments/' + expCaseDir)
	expProfile = {'motif': sys.argv[3],\
				'locus': 'loci/' + sys.argv[4] + '.bed',\
				'flankLength': int(sys.argv[5]),\
				'refGenomeDir': sys.argv[6]}
	motif = 
	locus = 'loci/' + sys.argv[4] + '.bed'
	flankLength = int(sys.argv[5])
	refGenomeDir = sys.argv[6]
	nCopyCount = int(sys.argv[7])
	nCopyList = []
	for i in range(8, 8 + nCopyCount):
		nCopyList.append(int(sys.argv[i]))

	with (expCaseDir + '/exp_profile.txt', 'w') as f:
		f.write('#### Experiment Profile ####')
		f.write('>>locus:\t', locus)
		f.write('>>motif:\t', motif)
		f.write('>>flankLength:\t', flankLength)
