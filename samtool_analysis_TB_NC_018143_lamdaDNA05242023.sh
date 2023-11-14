set -ue
sam_file=$1
#name=$2

filename_1=$(basename "$sam_file")

filename_1="${filename_1%%.sorted.bam}"

samtools sort -n $filename_1.sorted.bam  -o $filename_1.sorted_n.bam

samtools view -bf 0x2 $filename_1.sorted_n.bam | bedtools bamtobed -i stdin -bedpe > $filename_1\_PE.bedpe
cut -f 1,2,6,7,8,9 $filename_1\_PE.bedpe > $filename_1\_PE_2.bedpe
rm $filename_1\_PE.bedpe
bedtools getfasta -fi  ~/ref_genome/NC_018143_TB_lambda_genome.fna    -bed $filename_1\_PE_2.bedpe  -s  -name  > $filename_1\_PE_2_S.fa

bowtie2 -p 10   -f --very-sensitive-local    --ff -x ~/ref_genome/NC_018143_TB_lambda_genome  $filename_1\_PE_2_S.fa -S $filename_1\_PE_2_S.bowtie2mapping.sam

cat $filename_1\_PE_2_S.bowtie2mapping.sam  |cut -f 6 |sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1\_total_mapping.txt
samtools view -Sb $filename_1\_PE_2_S.bowtie2mapping.sam >$filename_1\_PE_2_S.bowtie2mapping.bam
samtools sort $filename_1\_PE_2_S.bowtie2mapping.bam > $filename_1\_PE_2_S.bowtie2mapping.sorted.bam
samtools index $filename_1\_PE_2_S.bowtie2mapping.sorted.bam 


bedtools intersect -abam  $filename_1\_PE_2_S.bowtie2mapping.sorted.bam -b  ~/ref_genome/TB_H37RV.RNA.bed  -v > $filename_1.rRNA_filtered.bam
samtools index $filename_1.rRNA_filtered.bam
samtools view  $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered.sam
cat $filename_1.rRNA_filtered.sam   |cut -f 6 | sed "s/M//g" |sed "s/\*$//g"| grep -vE "^-" | grep -vE "^0" >$filename_1.RNA_filtered_length.txt


rm $filename_1.rRNA_filtered.sam 

samtools sort $filename_1.rRNA_filtered.bam > $filename_1.rRNA_filtered_sorted.bam
samtools index $filename_1.rRNA_filtered_sorted.bam


samtools view -L ~/ref_genome/TB_H37RV.RNA.bed   $filename_1\_PE_2_S.bowtie2mapping.sorted.bam |cut -f 6 |sed "s/M//g" | grep -vE "^-" | grep -vE "^0" > $filename_1.rRNA_length.txt

file="txt"
if [ -d "$file" ]
then
	echo "$file exists and the new files will be added into txt file."
else
	mkdir txt
fi

mv $filename_1*.txt txt/


rm $filename_1\_PE_2_S.fa $filename_1\_PE_2_S.bowtie2mapping.sam



