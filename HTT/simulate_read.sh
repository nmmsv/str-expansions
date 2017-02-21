#!/bin/bash

mkdir -p fastq


for fasta_in in fasta/*
do
	tmp=${fasta_in%.fa}
	name=${tmp##*/}
	OUTDIR1=fastq/${name}.read1.fq
	OUTDIR2=fastq/${name}.read2.fq
	
	wgsim	$fasta_in \
		$OUTDIR1 \
		$OUTDIR2 \
		
done