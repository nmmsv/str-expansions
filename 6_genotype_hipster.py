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



expName = 'ATXN7_3'
locus = 'ATXN7'
nCopyList = [4, 5, 8, 10, 12, 14]
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
        --regions ' + repoDir + 'loci/' + locus + '.bed \
        --str-vcf ' + outDir + expName + '_' + nCopy + '.vcf.gz \
        --min-reads 5')

    os.system('gzip -d -f ' + outDir + expName + '_' + nCopy + '.vcf.gz')
