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



expName = 'ATXN7_err'
locus = 'ATXN7'
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
# ####
# expName = 'ATXN7_nCopy'                                                                                               
# nCopy = "list"                                                                                     
# numReads = '1000'         # number of read pairs                                                   
# error = '0'               # base error                                                             
# paramList = [0, 4, 8, 12, 16, 20, 24, 28, 40]               # nCopy list                           

# ####
expName = 'ATXN7_err'                                                                                               
nCopy = '8'                                                                                        
numReads = '1000'         # number of read pairs                                                   
error = "list"          # base error                                                               
paramList = [0, 0.0001, 0.001, 0.01, 0.03, 0.05, 0.08, 0.1] # error list                           

####
# expName = 'ATXN7_cov'                                                                                                 
# nCopy = '8'
# numReads = "list"       # number of read pairs                                                       
# error = '0'               # base error                                                               
# paramList = [30, 50, 100, 250, 500]                             # numRead (cov*10) list 

paramList = [str(i) for i in paramList]

if nCopy == "list" and numReads != "list" and error != "list":
        param = "nC"
elif nCopy != "list" and numReads == "list" and error != "list":
        param = "nR"
elif nCopy != "list" and numReads != "list" and error == "list":
        param = "er"
else:
        print 'Incorrect list assignment: (Only 1 list possible)'
        sys.exit()

print 'Creating new directories..'
outDir = repoDir + expName + '/Genotype_GATK/'
alignDir = repoDir + expName + '/aligned_read/'
TOOL = 'java -jar /storage/shahab/GenomeAnalysisTK.jar'
mkdir_p(outDir)
for pr in paramList:
    if param == "nC":
        nCopy = pr
    elif param == "nR":
        numRead = pr
    elif param == "er":
        error = pr
    else:
        print 'Wrong param = ', param
        sys.exit()
    alignCaseDir = alignDir + expName + '_' + param + '_' + pr + '/'
    outCaseDir = outDir + expName + '_' + param + '_' + pr + '/'
    mkdir_p(outCaseDir)
    print 'Calling GATK IndelRealigner for ', expName, ', nCopy=', nCopy + ', ' + param + ' ' + pr
    os.system(TOOL + ' ' \
        '-T IndelRealigner ' + \
        '-I ' + alignCaseDir + expName + '_' + nCopy + '.sorted.bam ' + \
        '-R ' + refGenomeDir + ' ' \
        '-targetIntervals ' + repoDir + 'loci' + '/target.list ' \
        '-o ' + outCaseDir + expName + '_' + nCopy + '_realigned.bam ')
    print 'Calling GATK HaplotypeCaller for ', expName, ', nCopy=', nCopy + ', ' + param + ' ' + pr
    os.system(TOOL + ' ' \
        '-T HaplotypeCaller ' + \
        '-I ' + outCaseDir + expName + '_' + nCopy + '_realigned.bam ' + \
        '-R ' + refGenomeDir + ' ' \
        '-L ' + repoDir + 'loci' + '/target.list ' \
#        '--minPruning 1 ' + \
#        '-minDanglingBranchLength 1 ' + \
#        '-maxAltAlleles 10 ' + \
#        '--maxNumHaplotypesInPopulation 1024 ' + \
#        '-ActProbThresh 0.0002 ' + \
#        '-indelHeterozygosity 0.5 ' + \
#        '--gcpHMM 100 ' + \
#        '-allowNonUniqueKmersInRef ' + \
#        '--kmerSize 35 ' + \
        '-o ' + outCaseDir + expName + '_' + nCopy + '.vcf ' + \
        '-bamout ' + outCaseDir + expName + '_' + nCopy + '.bamout.bam ')

