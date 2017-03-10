import subprocess
import errno    
import os
import sys
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


if len(sys.argv) == 1:
	print '### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads'
elif sys.argv[1] == "default":
	expName = 'HTT2'
	nCopyList = [21, 36]
	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
	repoDir = '/storage/nmmsv/str-expansions/'
	numThreads = 4
else:
	expName = sys.argv[1]
	refGenomeDir = sys.argv[2]
	repoDir = sys.argv[3]
	nCopyCount = int(sys.argv[4])
	nCopyList = []
	for i in range(5, 5 + nCopyCount):
		nCopyList.append(int(sys.argv[i]))
	numThreads = sys.argv[5 + nCopyCount]

print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList
print '## bwa mem parameters: ##'
print '# Number of threads: ', numThreads


print 

print '# Creating new directories..'
outDir = repoDir + expName + '/aligned_read/'
mkdir_p(outDir)

for nc in nCopyList:
	nCopy = str(nc)
	readCaseDir = repoDir + expName + '/simulated_read/' + expName + '_' + nCopy
	inFq1 = readCaseDir + '/' + expName + '_' + nCopy + '.read1.fq'
	inFq2 = readCaseDir + '/' + expName + '_' + nCopy + '.read2.fq'

	alignCaseDir = outDir + expName + '_' + nCopy + '/'
	mkdir_p(alignCaseDir)
	outFile = alignCaseDir + expName + '_' + nCopy

	#subprocess.call(['bwa', 'mem',\
	#		refGenomeDir, inFq1, inFq2,\
	#		'-R', '@RG\\tID:HTT\\tSM:16\\tLB:lb\\tPL:pl',\
	#		'>', outFile + '.sam'])
	print '# Executing bwa mem on ' + expName + ', nCopy=', nCopy	
	os.system('bwa mem ' + \
			refGenomeDir + ' ' +  inFq1 + ' ' + inFq2 +\
			' -R \'@RG\\tID:' + expName + '\\tSM:' + nCopy +\
			'\\tLB:lb\\tPL:pl\' > ' + outFile + '.sam')
	