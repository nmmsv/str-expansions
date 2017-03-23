import subprocess
import errno    
import os
import sys
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
# 	print '### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads'
# elif sys.argv[1] == "default":
# 	expName = 'HTT2'
# 	nCopyList = [21, 36]
# 	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
# 	repoDir = '/storage/nmmsv/str-expansions/'
# 	numThreads = 4
# else:
# 	expName = sys.argv[1]
# 	refGenomeDir = sys.argv[2]
# 	repoDir = sys.argv[3]
# 	nCopyCount = int(sys.argv[4])
# 	nCopyList = []
# 	for i in range(5, 5 + nCopyCount):
# 		nCopyList.append(int(sys.argv[i]))
# 	numThreads = sys.argv[5 + nCopyCount]

expName = sys.argv[1]
repoDir = sys.argv[2]

expCaseDir = repoDir + 'experiments/' + expName

try:
	with open(expCaseDir + '/profile.txt', 'r') as f:
		arg_dict = pickle.load(f)
except:
	raise

refGenomeDir = arg_dict['ref_gen_dir']
nCopyList = arg_dict['num_copy']
bamFilter = arg_dict['bam_filter']
numThreads = arg_dict['num_threads']
numReadsList = arg_dict['num_reads']
locus = 'loci/' + arg_dict['locus'] + '.bed'
l1 = arg_dict['read_len']

print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList
print '# Bam filter: ', bamFilter
print '# Number of read pairs: ', numReadsList
print '# Locus: ', locus
print '## bwa mem parameters: ##'
print '# Number of threads: ', numThreads


print 

print '# Creating new directories..'
expDir = repoDir + 'experiments/'
simReadDir = expDir + expName + '/simulated_read/'

outDir = expDir + expName + '/aligned_read/'
mkdir_p(outDir)

# First extract locus start and end
with open(repoDir + locus, 'r') as f:
	row = f.readline().split()
	chrom = row[0]
	locusStart = int(row[1])
	locusEnd = int(row[2])
# Setting the filter edges for left and right
# Looking for reads such that r1 <= leftFilter AND r2 >= rightFilter
leftFilter = locusStart - l1
rightFilter = locusEnd


for nc in nCopyList:
	nCopy = str(nc)
	readCaseDir = simReadDir + expName + '_nc' + nCopy
	alignCaseDir = outDir + expName + '_nc' + nCopy
	mkdir_p(alignCaseDir)
	for numReads in numReadsList:
		nRead = str(numReads)
		if len(numReadsList) > 1:
			subCaseName = 'nReads_' + nRead
			readSubCaseDir = readCaseDir + '/' + subCaseName
			alignSubCaseDir = alignCaseDir + '/' + subCaseName
			mkdir_p(alignSubCaseDir)
		else:
			readSubCaseDir = readCaseDir

		inFq1 = readSubCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead + '.read1.fq'
		inFq2 = readSubCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead + '.read2.fq'


		outFile = alignSubCaseDir + '/' + expName + '_nc' + nCopy + '_nr' + nRead

		#subprocess.call(['bwa', 'mem',\
		#		refGenomeDir, inFq1, inFq2,\
		#		'-R', '@RG\\tID:HTT\\tSM:16\\tLB:lb\\tPL:pl',\
		#		'>', outFile + '.sam'])
		print '# Executing bwa mem on ' + expName + ', nCopy =', nCopy, ' nReads =', nRead	
		os.system('bwa mem ' + \
				refGenomeDir + ' ' +  inFq1 + ' ' + inFq2 +\
				' -R \'@RG\\tID:' + expName + '\\tSM:' + nCopy +\
				'\\tLB:lb\\tPL:pl\' > ' + outFile + '.sam')
		os.system('samtools view -bT ' + refGenomeDir+\
				' ' + outFile + '.sam' + ' > '+\
				outFile + '.bam')
		os.system('samtools sort -o ' + outFile + '.sorted.bam '+\
				outFile + '.bam')
		os.system('samtools index ' + outFile + '.sorted.bam ' +\
				outFile + '.sorted.bai')
		os.system('samtools index ' + outFile + '.sorted.bam')

		if bamFilter:
			samDir = outFile + '.sam'
			newSamDir = outFile + '_flt.sam'
			tempDict = {}
			newSamFile = open(newSamDir, 'w')
			print 'Filtering ' + outFile + '.bam'
			with open(samDir, 'r') as samFile:
				for record in samFile:
					if record[0] == '@':
						newSamFile.write(record)
					else:
						row = record.split()
						if row[2] == chrom and int(row[3]) <= leftFilter and int(row[7]) >= rightFilter:
							newSamFile.write(record)

			newSamFile.close()

			os.system('samtools view -bT ' + refGenomeDir+\
					' ' + outFile + '_flt.sam' + ' > '+\
					outFile + '_flt.bam')
			os.system('samtools sort -o ' + outFile + '_flt.sorted.bam '+\
					outFile + '.bam')
			os.system('samtools index ' + outFile + '_flt.sorted.bam ' +\
					outFile + '_flt.sorted.bai')
		