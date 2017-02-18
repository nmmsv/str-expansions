#!/bin/bash
mkdir -p hpstr
if [ ! -f aligned/HTT_16.bam ]; then
	samtools view -bT /storage/resources/dbase/human/hs37d5/hs37d5.fa  aligned/HTT_16.sam > aligned/HTT_16.bam
fi

HipSTR --bams aligned/HTT_16.bam --fasta /storage/resources/dbase/human/hs37d5/hs37d5.fa  --regions HTT_STR_locus.bed --str-vcf hpstr/str_gts_16.vcf.gz
