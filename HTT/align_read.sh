#!/bin/bash

ref_genome=/storage/resources/dbase/human/hs37d5/hs37d5.fa
read1=fastq/HTT_str_16.read1.fq
read2=fastq/HTT_str_16.read2.fq
OUTDIR=aligned/HTT_16

bwa mem $ref_genome $read1 $read2 -R '@RG\tID:HTT\tSM:16\tLB:lb\tPL:pl' > $OUTDIR.sam

samtools view -bT /storage/resources/dbase/human/hs37d5/hs37d5.fa  $OUTDIR.sam > $OUTDIR.bam
samtools sort -o $OUTDIR.bam.sorted $OUTDIR.bam
samtools index $OUTDIR.bam.sorted $OUTDIR.bam.sorted.bai
#bwa mem -t 5 /storage/resources/dbase/human/hs37d5/hs37d5.fa fastq/HTT_str_16.read1.fq fastq/HTT_str_16/read2.fq
