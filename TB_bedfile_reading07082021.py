#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import re
from Bio.Seq import Seq
import glob, os
import statistics 
from statistics import mode 
import numpy

TSS_transcript_infor = {}
genome = {}
genome_seq_signal={}
cwd = os.getcwd()
input_file_name = os.path.splitext(str(sys.argv[1]))[0]
pattern = input_file_name.split(".bedpe")
input_file_name=pattern[0]









#home_directory = os.getenv("HOME")
for num in range (0,4471711):
    #genome[num] = [0,0,0,0,0,0,0,0]
    genome_seq_signal[num] = [num,0,0,0,0,0,0]
    
print ('finish the first step of data initialization\n')




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

exit()
    















exit()
genome_list={}
genome_infor={}
genome_input_file = open ("Bb_strain.fa", 'r')
for line in  genome_input_file:
    #print(line)
    if line == '\n':
        continue
    line = line.strip().split() #去除前后空格
    if re.search(r"^>", line[0]):
        
        
        #for num in range(1, (len(line)+1):
        new_line_2 = "_".join(str(i) for i in line)
        genome_name = new_line_2.replace(">", "")
        genome_name = genome_name.replace(",", "")
        gene_ID= line[0].replace(">", "")
        # (gene_ID, genome_name)
        genome_infor[gene_ID]=[gene_ID, genome_name]
    else:
        sequence = line[0]
        length = len(sequence)
        genome_infor[gene_ID].append(length)
        #print (gene_ID, length)
        genome_sub_name = gene_ID
        #print (genome_sub_name)
        vars()[genome_sub_name]={}
        for num in range (0,(len(sequence)+1000) ):
            vars()[genome_sub_name][int(num)]=[genome_sub_name, num,0,0,0,0,0,0,0] #0-3: TSS information, 4-6:positive strand: start, end, coverage, 
            #7-9: positive strand: end, start, coverage, #10, 11, 12 TTS information  #13, 14, 15, 16, 17, 18 genome annoation information
            #for num_add in range(0, sample_serial):
                #vars()[genome_sub_name][int(num)].append(int(0))
                #vars()[genome_sub_name][int(num)].append(int(0))
        


cwd = os.getcwd() 
input_file_name = os.path.splitext(str(sys.argv[1]))[0]
pattern = input_file_name.split(".bedpe")
input_file_name=pattern[0]

output_file=open(input_file_name+"_bed_reading.txt", 'w')
bedfile_full_directory = cwd+'/'+str(sys.argv[1])
print ("start to read the bedpe file: " + bedfile_full_directory + "\n")        
bed_input_file = open(bedfile_full_directory, 'r')


#bed_input_file= open ('B1_Borre_total_RD_1_S15_PE_2_S.bowtie2mapping.sorted.bedpe', 'r')
for line_bed in  bed_input_file:
    position_informaton_bed = line_bed.strip().split("\t")
        
    if (len(position_informaton_bed) >4) :
        (sub_genome_name, transcript_start, transcript_end, transcript_direction) = (position_informaton_bed[0], int(position_informaton_bed[1]), int(position_informaton_bed[2]), str(position_informaton_bed[5]))
        array_name= sub_genome_name+"_array"
        if transcript_direction == "+":
            
            genome_seq_signal[int(transcript_start)+1][3]+=1
            genome_seq_signal[int(transcript_end)+1][4]+=1
            for nucletide_position in range ((int(transcript_start)+1), (int(transcript_end)+1) ):
                genome_seq_signal[int(nucletide_position)][5]+=1
        if transcript_direction == "-":
            
            genome_seq_signal[int(transcript_start)][6]+=1
            genome_seq_signal[int(transcript_end)][7]+=1
            for nucletide_position in range ((int(transcript_start)), (int(transcript_end)) ):
                genome_seq_signal[int(nucletide_position)][8]+=1

bed_input_file.close



for sub_genome_name in genome_infor:
    print (sub_genome_name, genome_infor[sub_genome_name][2])
    if sub_genome_name == "Lambda_NEB": continue
    for nucletide in range (0, int(genome_infor[sub_genome_name][2])+1):
        new_line_2 = "\t".join(str(i) for i in genome_seq_signal[int(nucletide)])  + '\n' 
        output_file.write(new_line_2)
output_file.close()
