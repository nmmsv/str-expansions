exp=ATXN7_12_old_heatMap_d1000


python 0_create_profile.py 	--exp-name $exp \
							--locus ATXN7 \
							--motif GCA \
							--ref-allele-count 10 \
							--num-copy 3 7 10 16 20 30 40 55 70 85\
							--num-reads 150 360 600 900 1350 1800 2400 3000 \
							--ref-gen-dir /storage/resources/dbase/human/hs37d5/hs37d5.fa \
							--repo-dir /storage/nmmsv/str-expansions/ \
							--flank-len 3000 \
							--read-len 100 \
							--read-ins-mean 1000 \
							--read-ins-stddev 50 \
							--base-error 0 \
							--bam-filter True \
							--heat-map-limit 1



python 1_simulate_alt_genome.py $exp \
					/storage/nmmsv/str-expansions/


python 2.0.1_read_simulated_data_profile.py $exp \
					/storage/nmmsv/str-expansions/

python 3.0.1_align_read_profile.py $exp \
					/storage/nmmsv/str-expansions/

python vis_0.1_insert_size_profile.py $exp \
					/storage/nmmsv/str-expansions/ pdf