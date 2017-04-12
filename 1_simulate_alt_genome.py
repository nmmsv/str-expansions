import csv
import tempfile
import subprocess
import pickle
import errno    
import os
import sys

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
# print sys.argv


# if len(sys.argv) < 2:
# 	print '### Usage python 1_simulate_alt_genome.py expName locus motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ...'
# elif sys.argv[1] == "default":
# 	expName = 'HTT2'
# 	motif = 'CAG'
# 	locus = 'loci/HTT.bed'
# 	flankLength = 5000
# 	nCopyList = [21, 36]
# 	refGenomeDir = '/storage/resources/dbase/human/hs37d5/hs37d5.fa'
# 	repoDir = '/storage/nmmsv/str-expansions/'
# else:
# 	expName = sys.argv[1]
# 	motif = sys.argv[2]
# 	locus = 'loci/' + sys.argv[3] + '.bed'
# 	flankLength = int(sys.argv[4])
# 	refGenomeDir = sys.argv[5]
# 	repoDir = sys.argv[6]
# 	nCopyCount = int(sys.argv[7])
# 	nCopyList = []
# 	for i in range(8, 8 + nCopyCount):
# 		nCopyList.append(int(sys.argv[i]))

expName = sys.argv[1]
repoDir = sys.argv[2]

expCaseDir = repoDir + 'experiments/' + expName

try:
	with open(expCaseDir + '/profile.txt', 'r') as f:
		arg_dict = pickle.load(f)
except:
	raise
locus = 'loci/' + arg_dict['locus'] + '.bed'
motif = arg_dict['motif']
flankLength = arg_dict['flank_len']
refGenomeDir = arg_dict['ref_gen_dir']
nCopyList = arg_dict['num_copy']

print '## Simulating Alterred Genome ##'
print '# Experiment Name: ', expName
print '# Repeating Motif: ', motif
print '# Flank Length: ', flankLength
print '# Reference Genome Directory: ', refGenomeDir
print '# Repository Directory: ', repoDir
print '# List of nCopies for Generated Genomes: ', nCopyList

print '# Creating new directories..'
expDir = repoDir + 'experiments/'
mkdir_p(expDir + expName)
mkdir_p(expDir + expName + '/simulated_genome')
mkdir_p(expDir + expName + '/temp')
# STR_locus contains the information for the STR locus
# Format: Chr	Start	End	MotifLength	nCopies
with open (repoDir + locus, 'r') as f:
	for row in csv.reader(f, dialect='excel-tab'):
		print 'STR locus: ', row
		chrom = row[0]
		startLoc = int(row[1])
		endLoc = int(row[2])

beginFlank = startLoc - flankLength

endFlank = endLoc + flankLength

print '# Creating temp bed files..'
refBed = open(expDir + expName + '/temp/tempRefBed.bed', 'w')
refBed.write(chrom + '\t' + str(beginFlank) + '\t' + str(endFlank) + '\n')
refBed.close()

prefixBed = open(expDir + expName + '/temp/tempPrefixBed.bed', 'w')
prefixBed.write(chrom + '\t' + str(beginFlank) + '\t' + str(startLoc) + '\n')
prefixBed.close()

suffixBed = open(expDir + expName + '/temp/tempSuffixBed.bed', 'w')
suffixBed.write(chrom + '\t' + str(endLoc) + '\t' + str(endFlank) + '\n')
suffixBed.close()

refFastaDir = expDir + expName + '/temp/' + expName + '_ref.fa'
prefixFastaDir = expDir + expName + '/temp/' + expName + '_prefix.fa'
suffixFastaDir = expDir + expName + '/temp/' + expName + '_suffix.fa'
# Create blank fasta file if it doesn't exist:
open(refFastaDir, 'a').close()
open(prefixFastaDir, 'a').close()
open(suffixFastaDir, 'a').close()
# Extracting reference haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', expDir + expName + '/temp/tempRefBed.bed', \
	'-fo', refFastaDir])
with open (refFastaDir, 'r') as refFasta:
	refFasta.readline().strip()
	refHap = refFasta.readline().strip()

print '# Calling bedtools getfasta..'
# Generating prefix and suffix for altered (STR nCopy changed) haplotype
subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', expDir + expName + '/temp/tempPrefixBed.bed', \
	'-fo', prefixFastaDir])
with open (prefixFastaDir, 'r') as preFile:
	preFile.readline().strip()
	prefixHap = preFile.readline().strip()

subprocess.call(['bedtools', 'getfasta', \
	'-fi', refGenomeDir, \
	'-bed', expDir + expName + '/temp/tempSuffixBed.bed', \
	'-fo', suffixFastaDir])
with open (suffixFastaDir, 'r') as sufFile:
	sufFile.readline().strip()
	suffixHap = sufFile.readline().strip()

simGenDir = expDir + expName + '/simulated_genome/'

for nc in nCopyList:	
	with open(simGenDir + expName +'_' + str(nc) + '.fa', 'w') as f:
		f.write('>' + expName + '_' + str(nc) + '_altAllele\n')
		f.write(prefixHap + motif * nc + suffixHap + '\n')
		f.write('>' + expName + '_' + str(nc)+'_refAllele\n')
		f.write(refHap + '\n')
	print '# Created fasta file: ' + simGenDir + expName +'_' + str(nc) + '.fa'


