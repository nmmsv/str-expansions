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
outDir = repoDir + expName + '/Genotype_HipSTR/'
alignDir = repoDir + expName + '/aligned_read/'
mkdir_p(outDir)
for nc in nCopyList:
	nCopy = str(nc)
	print 'Calling HipSTR genotype tool for ', expName, ', nCopy=', nCopy
	os.system('HipSTR \
		--bams ' + alignDir + expName + '_' + nCopy + '/' + expName + '_' + nCopy + '.sorted.bam \
		--fasta ' + refGenomeDir + ' \
		--regions ' + repoDir + expName + '/' + expName + '_locus.bed \
		--str-vcf ' + outDir + expName + '_' + nCopy + '.vcf.gz')

