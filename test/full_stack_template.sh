

# expName=ATXN1_3
# motif=CAG
# disease=ATXN1
# flankLength=6000
# refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
# repoDir=/storage/nmmsv/str-expansions/
# dist=500
# stdDev=50
# numReads=100000

# expName=FMR1_4
# motif=GGC
# disease=FMR1
# flankLength=6000
# refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
# repoDir=/storage/nmmsv/str-expansions/
# l=120
# dist=1000
# stdDev=100
# numReads=100000


expName=ATXN7_5
motif=GCA
disease=ATXN7
flankLength=5000
refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
repoDir=/storage/nmmsv/str-expansions/
l=120
dist=500
stdDev=50
numReads=100000

cd /storage/nmmsv/str-expansions/

### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir 6 4 5 8 10 12 14

### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn dist stdDev l1 l2 numReads:
python 2_read_simulated_data.py $expName $motif $flankLength /storage/resources/dbase/human/hs37d5/hs37d5.fa /storage/nmmsv/str-expansions/ 6 4 5 8 10 12 14 $dist $stdDev $l $l $numReads

### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
python 3_align_read.py $expName $refGenomeDir $repoDir 6 4 5 8 10 12 14 4

python 5_filter_bam.py $expName $disease $refGenomeDir $repoDir $l 6 4 5 8 10 12 14

### Usage python vis_insert_size.py expName repoDir nCopyCount nCopy1 nCopy2 ... nCopyn extension
python vis_insert_size.py $expName $repoDir 6 4 5 8 10 12 14 .pdf


