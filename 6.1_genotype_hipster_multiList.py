import errno    
import os
import csv

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


expName = 'ATXN7_nCopy'
locus = 'ATXN7'
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'

#### 
nCopy = "list"
numReads = '1000'         # number of read pairs
error = '0'               # base error
paramList = [0, 4, 8, 12, 16, 20, 24, 28, 40]               # nCopy list

# ####
# nCopy = '8'
# numReads = '1000'         # number of read pairs
# error = "list"          # base error
# paramList = [0, 0.0001, 0.001, 0.01, 0.03, 0.05, 0.08, 0.09, 0.1, 0.11] # error list

# ####
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
outDir = repoDir + expName + '/Genotype_HipSTR/'
alignDir = repoDir + expName + '/aligned_read/'
mkdir_p(outDir)
for pr in paramList:
    if param == "nC":
        nCop = pr
        numR = numReads
        err  = error
    elif param == "nR":
        numR = pr
        nCop = nCopy
        err  = error
    elif param == "er":
        err  = pr
        numR = numReads
        nCop = nCopy
    else:
        print 'Wrong param = ', param
        sys.exit()
    print 'Calling HipSTR genotype tool for ', expName, ', nCopy=', nCop
    os.system('HipSTR \
        --bams ' + alignDir + expName + '_' + param + '_' + pr + '/' + expName + '_' + nCop + '.sorted.bam \
        --fasta ' + refGenomeDir + ' \
        --regions ' + repoDir + 'loci/' + locus + '.bed \
        --str-vcf ' + outDir + expName + '_' + param + '_' + pr + '.vcf.gz \
        --min-reads 2 --use-all-reads --use-unpaired')

    os.system('gzip -d -f ' + outDir + expName + '_' + param + '_' + pr + '.vcf.gz')


for pr in paramList:
    print '>>>>>>> param ', param, '\t', pr
    with open (outDir + expName + '_' + param + '_' + pr + '.vcf', 'r') as f:
        for line in csv.reader(f, dialect = 'excel-tab'):
            if len(line) > 5:
                print '\t'.join(line[0:5])
                print




