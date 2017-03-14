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



expName = 'ATXN7_4'
locus = 'ATXN7'
nCopyList = [5, 8, 10, 12, 14]
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
print 'Creating new directories..'
outDir = repoDir + expName + '/Genotype_GATK/'
alignDir = repoDir + expName + '/aligned_read/'
TOOL='java -jar /storage/shahab/GenomeAnalysisTK.jar'
mkdir_p(outDir)
for nc in nCopyList:
    nCopy = str(nc)
    print 'Calling GATK IndelRealigner for ', expName, ', nCopy=', nCopy
    os.system(TOOL + ' ' \
        '-T IndelRealigner ' + \
        '-I ' + alignDir + expName + '_' + nCopy + '/' + expName + '_' + nCopy + '.sorted.bam ' + \
        '-R ' + refGenomeDir + ' ' \
        '-targetIntervals ' + repoDir + expName + '/target.list ' \
        '-o ' + outDir + expName + '_' + nCopy + '_realigned.bam ')
    print 'Calling GATK HaplotypeCaller for ', expName, ', nCopy=', nCopy
    os.system(TOOL + ' ' \
        '-T HaplotypeCaller ' + \
        '-I ' + outDir + expName + '_' + nCopy + '_realigned.bam ' + \
        '-R ' + refGenomeDir + ' ' \
        '-L ' + repoDir + expName + '/target.list ' \
        '--minPruning 1 ' + \
        '-minDanglingBranchLength 1 ' + \
#        '-maxAltAlleles 10 ' + \
#        '--maxNumHaplotypesInPopulation 1024 ' + \
#        '-ActProbThresh 0.0002 ' + \
#        '-indelHeterozygosity 0.5 ' + \
#        '--gcpHMM 100 ' + \
        '-allowNonUniqueKmersInRef ' + \
        '--kmerSize 35 '
        '-o ' + outDir + expName + '_' + nCopy + '.vcf ' + \
        '-bamout ' + outDir + expName + '_' + nCopy + '.bamout.bam ')

