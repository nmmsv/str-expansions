import os, sys

if len(sys.argv) == 1:
	print '### Usage python vis_alignment.py expName repoDir refGenomeDir nCopy'
elif sys.argv[1] == "default":
	expName = 'HTT2'
	nCopy = str(21)
	repoDir = '/storage/nmmsv/str-expansions/'
	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
	ext = '.pdf'
else:
	expName = sys.argv[1]
	repoDir = sys.argv[2]
	refGenomeDir = sys.argv[3]
	nCopy = sys.argv[4]

print '## Aligning Reads from Simulated Data: ##'
print '# Experiment Name: ', expName
print '# Repository Directory: ', repoDir
print '# Regerence Genome Directory: ', refGenomeDir
print '# nCopy: ', nCopy

print 


alignDir = repoDir + expName + '/aligned_read/' + expName + '_' + nCopy + '/'
sortedBamDir = alignDir + expName + '_' + nCopy + '.sorted.bam'
os.system('samtools tview ' + sortedBamDir + ' ' + refGenomeDir)
