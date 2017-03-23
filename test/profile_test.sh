exp=ATXN7_12_heatMap_d1000


python 0_create_profile.py 	--exp-name $exp \
							--locus ATXN7 \
							--motif GCA \
							--ref-allele-count 10 \
							--num-copy 3 7 10 12 16 22 30 40 60 85\
							--num-reads 50 150 350 500 700 1000 2500 4000 8000 11000 \
							--ref-gen-dir /storage/resources/dbase/human/hs37d5/hs37d5.fa \
							--repo-dir /storage/nmmsv/str-expansions/ \
							--flank-len 4000 \
							--read-len 100 \
							--read-ins-mean 1000 \
							--read-ins-stddev 50 \
							--base-error 0 \
							--bam-filter True \
							--heat-map-limit 50



python 1_simulate_alt_genome.py $exp \
					/storage/nmmsv/str-expansions/


python 2.0.1_read_simulated_data_profile.py $exp \
					/storage/nmmsv/str-expansions/

python 3.0.1_align_read_profile.py $exp \
					/storage/nmmsv/str-expansions/

python vis_0.1_insert_size_profile.py $exp \
					/storage/nmmsv/str-expansions/ pdf