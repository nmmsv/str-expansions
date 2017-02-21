import os

expName = 'HTT2'
nCopy = str(36)
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
alignDir = repoDir + expName + '/aligned_read/' + expName + '_' + nCopy + '/'
sortedBamDir = alignDir + expName + '_' + nCopy + '.sorted.bam'
os.system('samtools tview ' + sortedBamDir + ' ' + refGenomeDir)