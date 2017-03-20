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


expName = 'ATXN7_err'
locus = ''
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'

#### 
# nCopy = "list"
# numReads = '1000'         # number of read pairs
# error = '0'               # base error
# paramList = [0, 4, 8, 12, 16, 20, 24, 28, 40]               # nCopy list

# ####
nCopy = '8'
numReads = '1000'         # number of read pairs
error = "list"          # base error
paramList = [0, 0.0001, 0.001, 0.01, 0.03, 0.05, 0.08, 0.09, 0.1, 0.11] # error list

# ####
# nCopy = '8'
# numReads = "list"       # number of read pairs
# error = '0'               # base error
# paramList = [30, 50, 100, 250, 500, 1000, 5000, 10000]                             # numRead (cov*10) list

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




expansionHunterPath = '/storage/nmmsv/ExpansionHunter/bin/ExpansionHunter'
print 'Creating new directories..'
outDir = repoDir + expName + '/Genotype_ExpansionHunter/'
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
    print 'Calling ExpansionHunter genotype tool for ', expName, ', nCopy=', nCop
    os.system(expansionHunterPath + ' \
        --bam ' + alignDir + expName + '_' + param + '_' + pr + '/' + expName + '_' + nCop + '.sorted.bam \
        --ref-fasta ' + refGenomeDir + ' \
        --repeat-specs ' + repoDir + 'loci/' + locus + ' \
        --vcf ' + outDir + expName + '_' + param + '_' + pr + '.vcf \
        --json ' + outDir + expName + '_' + param + '_' + pr + '.json \
        --log ' + outDir + expName + '_' + param + '_' + pr + '.log \
        --min-score 0.9 --min-baseq 20 --min-anchor-mapq 60 \
        --read-depth 1')

for pr in paramList:
    print '>>>>>>> param ', param, '\t', pr
    with open (outDir + expName + '_' + param + '_' + pr + '.vcf', 'r') as f:
        for line in csv.reader(f, dialect = 'excel-tab'):
            if len(line) > 5:
                print '\t'.join(line[0:5])
                print
