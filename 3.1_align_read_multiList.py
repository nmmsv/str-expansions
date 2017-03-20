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
	print '### Usage python 3_align_read.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN'
elif sys.argv[1] == "default":
	expName = 'HTT2'
	motif = 'CAG'
	flankLength = 5000
	nCopyList = [21, 36]
	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
	repoDir = '/storage/nmmsv/str-expansions/'
	dist = 500				# outer distance between the two pairs
	stdDev = 50				# Standard deviation
	l1 = 100				# length of first read
	l2 = 100				# length of second read
	numReads = 100000		# number of read pairs
else:
	expName = sys.argv[1]
	motif = sys.argv[2]
	flankLength = int(sys.argv[3])
	refGenomeDir = sys.argv[4]
	repoDir = sys.argv[5]
	nCopy = sys.argv[6]
	dist = sys.argv[7]			# outer distance between the two pairs
	stdDev = sys.argv[8]		# Standard deviation	
	l1 = sys.argv[9]			# length of first read
	l2 = sys.argv[9]			# length of second read
	numReads = sys.argv[10]		# number of read pairs
	error = sys.argv[11]		# base error

	lstCount = int(sys.argv[12])
	paramList = []
	for i in range(13, 13 + lstCount):
		paramList.append(sys.argv[i])

if nCopy == "list" and numReads != "list" and error != "list":
	param = "nC"
elif nCopy != "list" and numReads == "list" and error != "list":
	param = "nR"
elif nCopy != "list" and numReads != "list" and error == "list":
	param = "er"
else:
	print 'Incorrect list assignment: (Only 1 list possible)'
	sys.exit()

print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# nCopy: ', nCopy
print '# numReads: ', numReads
print '# error: ', error
print '## Parameter List: ', paramList
print '## bwa mem parameters: ##'
numThreads = 4
print '# Number of threads: Fixed in this version', numThreads


print 

print '# Creating new directories..'
expDir = repoDir + 'experiments/'
outDir = expDir + expName + '/aligned_read/'
mkdir_p(outDir)
simReadDir = expDir + expName + '/simulated_read/'
for pr in paramList:

	if param == "nC":
		nCopy = pr
	elif param == "nR":
		numRead = pr
	elif param == "er":
		error = pr
	else:
		print 'Wrong param = ', param
		sys.exit()

	readCaseDir = simReadDir + expName + '_' + param + '_' + pr

	inFq1 = readCaseDir + '/' + expName + '_' + nCopy + '.read1.fq'
	inFq2 = readCaseDir + '/' + expName + '_' + nCopy + '.read2.fq'

	alignCaseDir = outDir + expName + '_' + param + '_' + pr + '/'
	mkdir_p(alignCaseDir)
	outFile = alignCaseDir + expName + '_' + nCopy

	#subprocess.call(['bwa', 'mem',\
	#		refGenomeDir, inFq1, inFq2,\
	#		'-R', '@RG\\tID:HTT\\tSM:16\\tLB:lb\\tPL:pl',\
	#		'>', outFile + '.sam'])
	print '# Executing bwa mem on ' + expName + ', nCopy=', nCopy + ', ' + param + ' ' + pr 
	os.system('bwa mem ' + \
			refGenomeDir + ' ' +  inFq1 + ' ' + inFq2 +\
			' -t ' + str(numThreads) + ' -R \'@RG\\tID:' + expName + '\\tSM:' + nCopy +\
			'\\tLB:lb\\tPL:pl\' > ' + outFile + '.sam')
	os.system('samtools view -bT ' + refGenomeDir+\
			' ' + outFile + '.sam' + ' > '+\
			outFile + '.bam')
	os.system('samtools sort -o ' + outFile + '.sorted.bam '+\
			outFile + '.bam')
	os.system('samtools index ' + outFile + '.sorted.bam ' +\
			outFile + '.sorted.bai')
	os.system('samtools index ' + outFile + '.sorted.bam')
	