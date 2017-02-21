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
nCopyList = [36]
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
print 'Creating new directories..'
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
	print 'Executing bwa mem on ' + expName + ', nCopy=', nCopy	
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