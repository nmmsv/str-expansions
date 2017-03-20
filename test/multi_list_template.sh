# expName=ATXN7_err
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
# nCopyList="1 8"	# for 1_simulate_alt_genome.py
# nCopy=8
# lst="10 0 0.0001 0.001 0.01 0.03 0.05 0.08 0.09 0.1 0.11"

# cd /storage/nmmsv/str-expansions/

# ### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
# python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

# ### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
# python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

# ### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
# python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst


###########
expName=ATXN7_nCopy3
motif=GCA
disease=ATXN7
flankLength=1000
refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
repoDir=/storage/nmmsv/str-expansions/
l=100
dist=500
stdDev=50
numReads=5000
error=0
nCopyList="9 0 4 8 12 16 20 24 28 40"	# for 1_simulate_alt_genome.py
nCopy="list"
lst="9 0 4 8 12 16 20 24 28 40"

cd /storage/nmmsv/str-expansions/

### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

##########

# expName=ATXN7_cov
# motif=GCA
# disease=ATXN7
# flankLength=1000
# refGenomeDir=/storage/resources/dbase/human/hs37d5/hs37d5.fa
# repoDir=/storage/nmmsv/str-expansions/
# l=100
# dist=500
# stdDev=50
# numReads="list"
# error=0
# nCopyList="1 8"	# for 1_simulate_alt_genome.py
# nCopy=8
# lst="8 30 50 100 250 500 1000 5000 10000"

# cd /storage/nmmsv/str-expansions/

# ### Usage python 1_simulate_alt_genome.py expName motif flankLength refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn: 
# python 1_simulate_alt_genome.py $expName $motif $disease $flankLength $refGenomeDir $repoDir $nCopyList

# ### Usage python 2_read_simulated_data.py expName motif flankLength refGenomeDir repoDir nCopy dist stdDev l1 numReads error lstCount lst1 lst2 ... lstN
# python 2.1_read_simulated_data_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst

# ### Usage python 3_align_read.py expName refGenomeDir repoDir nCopyCount nCopy1 nCopy2 ... nCopyn numThreads
# python 3.1_align_read_multiList.py $expName $motif $flankLength $refGenomeDir $repoDir $nCopy $dist $stdDev $l $numReads $error $lst