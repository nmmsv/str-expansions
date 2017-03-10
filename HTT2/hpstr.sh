#!/bin/bash


HipSTR --bams aligned_read/HTT2_21/HTT2_21.sorted.bam --fasta /storage/resources/dbase/human/hs37d5/hs37d5.fa  --regions HTT2_locus.bed --str-vcf Genotype_HipSTR/str_gts_16.vcf.gz
