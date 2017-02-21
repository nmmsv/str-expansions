import csv
import tempfile
import subprocess

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
motif = 'CAG'
flankLength = 5000
nCopyList = [21, 36]
refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
repoDir = '/storage/nmmsv/str-expansions/'
print 'Creating new directories..'
mkdir_p(repoDir + expName)
mkdir_p(repoDir + expName + '/simulated_genome')
mkdir_p(repoDir + expName + '/temp')
# STR_locus contains the information for the STR locus
# Format: Chr	Start	End	MotifLength	nCopies
with open (repoDir + expName + '/' + expName + '_locus.bed', 'r') as f:
	for row in csv.reader(f, dialect='excel-tab'):
		print 'STR locus: ', row
		chrom = int(row[0])
		startLoc = int(row[1])
		endLoc = int(row[2])

beginFlank = startLoc - flankLength

endFlank = endLoc + flankLength

print 'Creating temp bed files..'
refBed = open(repoDir + expName + '/temp/tempRefBed.bed', 'w')
refBed.write(str(chrom) + '\t' + str(beginFlank) + '\t' + str(endFlank) + '\n')
refBed.close()

prefixBed = open(repoDir + expName + '/temp/tempPrefixBed.bed', 'w')
prefixBed.write(str(chrom) + '\t' + str(beginFlank) + '\t' + str(startLoc) + '\n')
prefixBed.close()

suffixBed = open(repoDir + expName + '/temp/tempSuffixBed.bed', 'w')
suffixBed.write(str(chrom) + '\t' + str(endLoc) + '\t' + str(endFlank) + '\n')
suffixBed.close()

refFastaDir = repoDir + expName + '/temp/' + expName + '_ref.fa'
prefixFastaDir = repoDir + expName + '/temp/' + expName + '_prefix.fa'
suffixFastaDir = repoDir + expName + '/temp/' + expName + '_suffix.fa'
# Create blank fasta file if it doesn't exist:
open(refFastaDir, 'a').close()
open(prefixFastaDir, 'a').close()
open(suffixFastaDir, 'a').close()
# Extracting reference haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', repoDir + expName + '/temp/tempRefBed.bed', \
	'-fo', refFastaDir])
with open (refFastaDir, 'r') as refFasta:
	refFasta.readline().strip()
	refHap = refFasta.readline().strip()

print 'Calling bedtools getfasta..'e
# Generating prefix and suffix for altered (STR nCopy changed) haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', repoDir + expName + '/temp/tempPrefixBed.bed', \
	'-fo', prefixFastaDir])
with open (prefixFastaDir, 'r') as preFile:
	preFile.readline().strip()
	prefixHap = preFile.readline().strip()

subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', repoDir + expName + '/temp/tempSuffixBed.bed', \
	'-fo', suffixFastaDir])
with open (suffixFastaDir, 'r') as sufFile:
	sufFile.readline().strip()
	suffixHap = sufFile.readline().strip()

simGenDir = repoDir + expName + '/simulated_genome/'

for nc in nCopyList:	
	with open(simGenDir + expName +'_' + str(nc) + '.fa', 'w') as f:
		f.write('>' + expName + '_' + str(nc) + '_altAllele\n')
		f.write(prefixHap + motif * nc + suffixHap + '\n')
		f.write('>' + expName + '_' + str(nc)+'_refAllele\n')
		f.write(refHap + '\n')
	print 'Created fasta file: ' + simGenDir + expName +'_' + str(nc) + '.fa'


