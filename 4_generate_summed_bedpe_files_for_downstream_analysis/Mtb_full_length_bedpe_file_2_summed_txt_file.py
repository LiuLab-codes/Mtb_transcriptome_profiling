#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path


#The input file should be full-length bedpe files, and the output will be a summed text file. The output will detail the following information from the first column onwards: 'genome_position', 'positive_strand_5'end_intensity', 'positive_strand_3'end_intensity', 'positive_strand_RNA_coverage', 'negative_strand_3'end_intensity', 'negative_strand_5'end_intensity', and 'negative_strand_RNA_coverage'.




genome_seq_signal={}
cwd = os.getcwd()

for num in range (0,4471711):

    genome_seq_signal[num] = [num,0,0,0,0,0,0]
 

input_file_name = os.path.splitext(str(sys.argv[1]))[0]
bedfile_full_directory = cwd+'/'+str(sys.argv[1])
pattern = input_file_name.split("_PE")
input_file_name_s=pattern[0]
print ("start to read the bedpe file: " + bedfile_full_directory + "\n")        
bed_input_file = open(bedfile_full_directory, 'r')
output_file=open(input_file_name_s+"_bed_reading.txt", 'w')
for line_bed in  bed_input_file:
    position_informaton_bed = line_bed.strip().split("\t")
        
    if (len(position_informaton_bed) >4) :
        (sub_genome_name, transcript_start, transcript_end, transcript_direction) = (position_informaton_bed[0], int(position_informaton_bed[1]), int(position_informaton_bed[2]), str(position_informaton_bed[5]))
        #array_name= sub_genome_name+"_array"
        if transcript_direction == "+":
            
            genome_seq_signal[int(transcript_start)+1][1]+=1
            genome_seq_signal[int(transcript_end)][2]+=1
            for nucletide_position in range ((int(transcript_start)+1), (int(transcript_end)+1) ):
                genome_seq_signal[int(nucletide_position)][3]+=1
        if transcript_direction == "-":
            
            genome_seq_signal[int(transcript_start)+1][4]+=1
            genome_seq_signal[int(transcript_end)][5]+=1
            for nucletide_position in range ((int(transcript_start)+1), (int(transcript_end)+1) ):
                genome_seq_signal[int(nucletide_position)][6]+=1

bed_input_file.close



for position in range (0,4471711):
    new_line_2 = "\t".join(str(i) for i in genome_seq_signal[int(position)])  + '\n' 
    output_file.write(new_line_2)
output_file.close()


    