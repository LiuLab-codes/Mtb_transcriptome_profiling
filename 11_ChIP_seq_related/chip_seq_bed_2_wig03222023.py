#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import re
from Bio.Seq import Seq
import glob, os


# The input files for this script should be the *.sam files generated from the previous ChIP-seq genome mapping step. The output file *.wig can be visualized using IGV with the same reference genome. Additionally, the *bed_reading.txt file is used for downstream analysis (column 1: genome position, column2: intensity of ChIP-seq signal at each position).

genome_seq_signal={}

genome = {}
for num in range (0,4471711):

    genome_seq_signal[num] = [num,0]
    
print ('finish the first step of data initialization\n')

input_file= str(sys.argv[1])
print(input_file)
input_file_name = os.path.splitext(str(sys.argv[1]))[0]
input_file_name_n= input_file_name+ "_n.bam"

os.system("samtools sort -n " + input_file + " -o " + input_file_name_n )
bedpe_name=input_file_name+"_PE.bedpe"
os.system("samtools view -bf 0x2 " + input_file_name_n + " | bedtools bamtobed -i stdin -bedpe > " + bedpe_name)
os.system("rm " + input_file_name_n)
print(input_file_name)


bedpe_input= open(bedpe_name, 'r')

for line in  bedpe_input:
    position_informaton = line.strip().split("\t")
    if len(position_informaton) >6:
        for nucletide_position in range (int(position_informaton[1]), int(position_informaton[5])):
            genome_seq_signal[int(nucletide_position)][1]+=1


bedpe_input.close()
            
        
output_file_wig= open(input_file_name+ '_Mtb_chip_seq.wig',  'w')
output_file_txt= open(input_file_name+ '_Mtb_chip_seq_bed_reading.txt',  'w')
output_file_wig.write("track name=\"" + input_file_name +"\" " +"color=98,0,234 altColor=0,0,255 graphType=bar viewLimits=0:100 \nfixedStep chrom=gi|561108321|ref|NC_018143.2| start=1 step=1\n")

for position in genome_seq_signal:
    output_file_wig.write(str(genome_seq_signal[position][1])+"\n")
    output_file_txt.write(str(position) +"\t" + str(genome_seq_signal[position][1])+"\n")
    #output_file_n.write(str(genome[position][3][6])+"\n")
    
output_file_wig.close()
output_file_txt.close()
os.system("rm " + bedpe_name)

        
    
