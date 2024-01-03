set -ue

#“This script can automatically detect the paired *fastq.gz files in the same directory. However, you need to specify the reference genome at line 13. To use the script, run the following command: sh Mtb_ChIP_data_bowtie2_mapping_pipline.sh”

for r1 in *_R1_001.fastq.gz 
do
	
r2=${r1/_R1/_R2} 
out=${r1%%_R1_001.fastq.gz}.sam 
bam_out=${r1%%_R1_001.fastq.gz}.bam 
bam_out_sort=${r1%%_R1_001.fastq.gz}.sorted.bam 

(bowtie2 -p4 --no-mixed  -x ~/ref_genome/TB_H37RV  -1 $r1 -2 $r2 -S $out) 2> $out.log     #The reference genome directory is subjected to change. 
samtools view -u  $out | samtools sort -o $bam_out_sort

samtools index $bam_out_sort

rm $out 
done
