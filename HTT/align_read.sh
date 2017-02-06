#!/bin/bash

ref_genome=/storage/resources/dbase/human/hs37d5/hs37d5.fa
read1=fastq/HTT_str_16.read1.fq
read2=fastq/HTT_str_16.read2.fq
OUTDIR=out.asm

bwa mem $ref_genome $read1 $read2 
#bwa mem -t 5 /storage/resources/dbase/human/hs37d5/hs37d5.fa fastq/HTT_str_16.read1.fq fastq/HTT_str_16/read2.fq
