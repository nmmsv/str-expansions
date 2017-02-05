#!/bin/bash
# Extracting target STR loci from reference genome.

target_loci=$1
output_fasta=$2
ref_genome=/storage/resources/dbase/human/hs37d5/hs37d5.fa

bedtools getfasta -fi $ref_genome -bed $target_loci -fo $output_fasta

