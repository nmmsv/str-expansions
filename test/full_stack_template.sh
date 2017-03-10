# First Create following file (and directory):
# 	repoDir/expName/expName_locus.bed
expName=ATXN1_2
motif=CAG
disease=ATXN1
flankLength=5000
refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
repoDir=/storage/nmmsv/str-expansions/
dist=1500
stdDev=100
numReads=100000

### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir 6 10 20 40 80 120 150

### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn dist stdDev l1 l2 numReads:
python 2_read_simulated_data.py $expName $motif $flankLength /storage/resources/dbase/human/hs37d5/hs37d5.fa /storage/nmmsv/str-expansions/ 6 10 20 40 80 120 150 $dist $stdDev 100 100 $numReads

### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
python 3_align_read.py $expName /storage/resources/dbase/human/hs37d5/hs37d5.fa /storage/nmmsv/str-expansions/ 6 10 20 40 80 120 150 4

### Usage python vis_insert_size.py expName repoDir nCopyCount nCopy1 nCopy2 ... nCopyn extension
python vis_insert_size.py $expName /storage/nmmsv/str-expansions/ 6 10 20 40 80 120 150 .pdf


