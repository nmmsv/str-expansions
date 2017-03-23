import subprocess
import sys

import errno    
import os
import pickle

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


# if len(sys.argv) == 1:
# 	print '### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn dist stdDev l1 l2 numReads'
# elif sys.argv[1] == "default":
# 	expName = 'HTT2'
# 	motif = 'CAG'
# 	flankLength = 5000
# 	nCopyList = [21, 36]
# 	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
# 	repoDir = '/storage/nmmsv/str-expansions/'
# 	dist = 500				# outer distance between the two pairs
# 	stdDev = 50				# Standard deviation
# 	l1 = 100				# length of first read
# 	l2 = 100				# length of second read
# 	numReads = 100000		# number of read pairs
# else:
# 	expName = sys.argv[1]
# 	motif = sys.argv[2]
# 	flankLength = int(sys.argv[3])
# 	refGenomeDir = sys.argv[4]
# 	repoDir = sys.argv[5]
# 	nCopyCount = int(sys.argv[6])
# 	nCopyList = []
# 	for i in range(7, 7 + nCopyCount):
# 		nCopyList.append(int(sys.argv[i]))
# 	dist = sys.argv[7 + nCopyCount]			# outer distance between the two pairs
# 	stdDev = sys.argv[8 + nCopyCount]		# Standard deviation
# 	l1 = sys.argv[9 + nCopyCount]			# length of first read
# 	l2 = sys.argv[10 + nCopyCount]			# length of second read
# 	numReads = sys.argv[11 + nCopyCount]	# number of read pairs



expName = sys.argv[1]
repoDir = sys.argv[2]

expCaseDir = repoDir + 'experiments/' + expName

try:
	with open(expCaseDir + '/profile.txt', 'r') as f:
		arg_dict = pickle.load(f)
except:
	raise
locus = 'loci/' + arg_dict['locus'] + '.bed'
motif = arg_dict['motif']
flankLength = arg_dict['flank_len']
refGenomeDir = arg_dict['ref_gen_dir']
nCopyList = arg_dict['num_copy']
dist = arg_dict['read_ins_mean']
stdDev = arg_dict['read_ins_stddev']
l1 = arg_dict['read_len']
l2 = l1
baseError = arg_dict['base_error']
numReadsList = arg_dict['num_reads']


print '## Generating Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Repeating Motif: ', motif
print '# Flank Length: ', flankLength
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList
print '## wgsim parameters: ##'
print '# Mean outer distance between two pairs (-d): ', dist
print '# Standard Deviation of outer distance (-s): ', stdDev
print '# Length of first read (-1): ', l1
print '# Length of second read (-2): ', l2
print '# Number of read pairs (-N): ', numReadsList
expDir = repoDir + 'experiments/'
simReadDir = expDir + expName + '/simulated_read/'
print '# Creating new directories..'
mkdir_p(simReadDir)

for nc in nCopyList:
	nCopy = str(nc)
	caseDir = simReadDir + expName + '_nc' + nCopy
	mkdir_p(caseDir)
	for numReads in numReadsList:
		nRead = str(numReads)
		if len(numReadsList) > 1:
			subCaseName = 'nReads_' + nRead
			subCaseDir = caseDir + '/' + subCaseName
			mkdir_p(subCaseDir)
		else:
			subCaseDir = caseDir
		
		inFa = expDir + expName + '/simulated_genome/' + expName + '_' + nCopy + '.fa'
		outFq1 = subCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead + '.read1.fq'
		outFq2 = subCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead + '.read2.fq'
		print '# Generating files: ', outFq1
		subprocess.call(['wgsim',\
				'-d', str(dist),\
				'-s', str(stdDev),\
	            '-e', str(baseError),\
	#                        '-r', str(0.0),\
	                        '-1', str(l1),\
				'-2', str(l2),\
				'-N', str(numReads),\
				inFa, outFq1, outFq2])

