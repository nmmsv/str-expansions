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



expName = 'ATXN7_5'
locus = 'EH_ATXN7.json'
locus = ''
nCopyList = [4, 5, 8, 10, 12, 14]
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
expansionHunterPath = '/storage/nmmsv/ExpansionHunter/bin/ExpansionHunter'
print 'Creating new directories..'
outDir = repoDir + expName + '/Genotype_ExpansionHunter/'
alignDir = repoDir + expName + '/aligned_read/'
mkdir_p(outDir)
for nc in nCopyList:
    nCopy = str(nc)
    print 'Calling ExpansionHunter genotype tool for ', expName, ', nCopy=', nCopy
    os.system(expansionHunterPath + ' \
        --bam ' + alignDir + expName + '_' + nCopy + '/' + expName + '_' + nCopy + '.sorted.bam \
        --ref-fasta ' + refGenomeDir + ' \
        --repeat-specs ' + repoDir + 'loci/' + locus + ' \
        --vcf ' + outDir + expName + '_' + nCopy + '.vcf \
        --json ' + outDir + expName + '_' + nCopy + '.json \
        --log ' + outDir + expName + '_' + nCopy + '.log \
        --min-score 0.1 --min-baseq 1 --min-anchor-mapq 1 \
        --read-depth 1')

