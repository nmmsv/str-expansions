import subprocess



import errno    
import os

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

expName = 'HTT2'
motif = 'CAG'
flankLength = 5000
nCopyList = [21, 36]
dist = 500				# outer distance between the two pairs
stdDev = 50				# Standard deviation
l1 = 100				# length of first read
l2 = 100				# length of second read
numReads = 100000		# number of read pairs

repoDir = '/storage/nmmsv/str-expansions/'
simReadDir = repoDir + expName + '/simulated_read/'
print 'Creating new directories..'
mkdir_p(simReadDir)

for nc in nCopyList:
	nCopy = str(nc)
	caseDir = simReadDir + expName + '_' + nCopy
	inFa = repoDir + expName + '/simulated_genome/' + expName + '_' + nCopy + '.fa'
	outFq1 = caseDir + '/' + expName + '_' + nCopy + '.read1.fq'
	outFq2 = caseDir + '/' + expName + '_' + nCopy + '.read2.fq'
	mkdir_p(caseDir)

	subprocess.call(['wgsim',\
			'-d', str(dist),\
			'-s', str(stdDev),\
			'-1', str(l1),\
			'-2', str(l2),\
			'-N', str(numReads),\
			inFa, outFq1, outFq2])
	with open(caseDir + '/read_info.txt', 'w') as f:
		f.write('Simulated read using \'wgsim\' for experiment ' + expName + ', nCopy = ' + nCopy + '\n')
		f.write('Parameters:\n')
		f.write('-d\t' + str(dist) + '\t\tOuter distance between two pairs\n')
		f.write('-s\t' + str(stdDev) + '\t\tStandard deviation\n')
		f.write('-N\t' + str(numReads) + '\tNumber of reads\n')
		f.write('-1\t' + str(l1) + '\t\tLength of first Read\n')
		f.write('-2\t' + str(l2) + '\t\tLength of second Read\n')

