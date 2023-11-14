#!/bin/bash

for r1 in *_split_R1.fa 
do
	
r2=${r1/_R1/_R2} 
out=${r1%%_split_R1.fa}.sam 
bam_out=${r1%%_split_R1.fa}.bam 
bam_out_sort=${r1%%_split_R1.fa}.sorted.bam 
timestamp=$(date +%s)
(bowtie2 -p 8   -N 1 --no-mixed   -f --very-sensitive-local   -X 10000 --ff -x ~/ref_genome/NC_018143_TB_lambda_genome  -1 $r1 -2 $r2 -S $out) 2> $out\_$timestamp.log
samtools view -u  $out | samtools sort -o $bam_out_sort
# -N 1 --no-mixed  -f --very-sensitive-local   -X 10000 --ff
samtools index $bam_out_sort

rm $out 
done
