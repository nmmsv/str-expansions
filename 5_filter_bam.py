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
	print '### Usage python 5_filter_bam.py expName locus refGenomeDir repoDir l1 nCopyCount nCopy1 nCopy2 ... nCopyn'
	sys.exit()
elif sys.argv[1] == "default":
	expName = 'HTT2'
	locus = 'loci/HTT.bed'
	nCopyList = [21, 36]
	repoDir = '/storage/nmmsv/str-expansions/'
	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
	l1 = 100				# length of first read
else:
	expName = sys.argv[1]
	locus = 'loci/' + sys.argv[2] + '.bed'
	refGenomeDir = sys.argv[3]
	repoDir = sys.argv[4]
	l1 = int(sys.argv[5])
	nCopyCount = int(sys.argv[6])
	nCopyList = []
	for i in range(7, 7 + nCopyCount):
		nCopyList.append(int(sys.argv[i]))

print '## Filtering bam file: Only keeping reads spanning the STR locus ##'
print '# Experiment Name: ', expName
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# Length of read1: ', l1
print '# List of nCopies for Generated Genomes: ', nCopyList


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
	outFile = repoDir + expName + '/aligned_read/' + expName + '_' + nCopy + '/' + expName + '_' + nCopy
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