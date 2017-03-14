# expName=ATXN7_multi_test
# motif=GCA
# disease=ATXN7
# flankLength=1000
# refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
# repoDir=/storage/nmmsv/str-expansions/
# l=100
# dist=500
# stdDev=50
# numReads=1000
# error="list"
# nCopyList="1 5"	# for 1_simulate_alt_genome.py
# nCopy=5
# lst="5 0 0.0001 0.001 0.01 0.1"

# cd /storage/nmmsv/str-expansions/

# ### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
# python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

# ### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
# python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

# ### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
# python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst


###########
# expName=ATXN7_multi_test_2
# motif=GCA
# disease=ATXN7
# flankLength=1000
# refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
# repoDir=/storage/nmmsv/str-expansions/
# l=100
# dist=500
# stdDev=50
# numReads=1000
# error=0
# nCopyList="3 4 8 16"	# for 1_simulate_alt_genome.py
# nCopy="list"
# lst="3 4 8 16"

# cd /storage/nmmsv/str-expansions/

# ### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
# python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

# ### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
# python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

# ### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
# python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

##########

expName=ATXN7_multi_test_3
motif=GCA
disease=ATXN7
flankLength=1000
refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
repoDir=/storage/nmmsv/str-expansions/
l=100
dist=500
stdDev=50
numReads="list"
error=0
nCopyList="1 5"	# for 1_simulate_alt_genome.py
nCopy=5
lst="3 100 1000 10000"

cd /storage/nmmsv/str-expansions/

### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst