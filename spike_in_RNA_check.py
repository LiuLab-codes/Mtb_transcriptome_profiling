#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import numpy
import scipy
import scipy.stats
import scipy.optimize
import math
import re
import logging
#from rdp import rdp
import matplotlib.pyplot as plt
from numpy import cov
import time

TSS_transcript_infor = {}
genome = {}
genome_seq_signal={}
cwd = os.getcwd()
#input_file_name = os.path.splitext(str(sys.argv[1]))[0]
home_directory = os.getenv("HOME")
for num in range (0,4471711):
#for num in range (0,271711):
    genome[num] = ["",0,0,"","","",[],[]]
    #genome_seq_signal[num] = [num,0,0,0,0,0,0]
    #if num > 220000 : break
    
print ('finish the first step of data initialization\n')

genome_annotation_input_file = open(home_directory + '/ref_genome/TB_H37RV.ref_genome_annoation06082020.txt', 'r')
for line in  genome_annotation_input_file:
    
    position_informaton = line.strip().split()
    #if int(position_informaton[0]) > 210000 : break
        #break
	
	#7	+	RVBD_0001:+:1:1524	0	null:+:0:0:7	RVBD_0002:+:2052:3260:2045
    #new_line = str(line) + '\n'
    #output_file.write(new_line)
    #length_infor = len(position_informaton)
    #print (length_infor )
    if (len(position_informaton) >5) :
        genome[int(position_informaton[0])][0] = [position_informaton[0], position_informaton[1], position_informaton[2], position_informaton[3], position_informaton[4], position_informaton[5]]
        #print(genome[int(position_informaton[0])][0])

genome_annotation_input_file.close()


print ("finish the step of loading the genome annotation information\n")


TSS_input_file =  open('TSS_collection_verified_by_total_RNA07082021.txt', 'r')

for line in  TSS_input_file:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: continue
        #if int(position_informaton[0]) >21000:
            #continue
    
        if (len(position_informaton) >3):
            if (str(position_informaton[1]) == "+"):
                genome[int(position_informaton[0])][1]=int(position_informaton[2])
                
                
            elif (str(position_informaton[1]) == "-"):
                genome[int(position_informaton[0])][2]=int(position_informaton[2])
        
    #print (position_informaton[0], position_informaton[1], position_informaton[2])
TSS_input_file.close()


print ("finish the step of loading the TSS information\n")

input_file_name = str(sys.argv[1])
pattern = input_file_name.split("_bed_reading")
input_file_name_s=pattern[0]
#print(input_file_name_s)
out_put_file= open('spike_in_RNA_check.txt', 'a')
reads_number_total, spike_in_RNA_reads_number=(0,0)
spike_in_RNA_coverage=1
RNA_coverage=0
spike_in_RNA_s=0
bed_read_file_1 = open(input_file_name, 'r')
for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        
        if   1471585 < int(position_informaton[0]) < 1477050 : 
            genome[int(position_informaton[0])][3]=[int(position_informaton[0]), 0, 0,0,0,0,0]
        
        else:
            genome[int(position_informaton[0])][3] = position_informaton
        
        
        
            if int(position_informaton[0]) <=4411709: RNA_coverage += int(genome[int(position_informaton[0])][3][3])
            if int(position_informaton[0]) <=4411709:
                reads_number_total+=1 
                RNA_coverage += int(genome[int(position_informaton[0])][3][6])
            if int(position_informaton[0]) > 4411709:
                spike_in_RNA_reads_number+=1
                spike_in_RNA_coverage +=int(genome[int(position_informaton[0])][3][3])
            if   2500386 < int(position_informaton[0]) < 2500832 :
                print(position_informaton[6])
                spike_in_RNA_s+=int(genome[int(position_informaton[0])][3][6])
                
                 
        
        

        
bed_read_file_1.close()   
spike_in_RNA_coverage=round(spike_in_RNA_coverage/60000,2 )+1
RNA_coverage_normalized= round(RNA_coverage/4456782, 2)

print(input_file_name_s,RNA_coverage_normalized, spike_in_RNA_coverage, reads_number_total,spike_in_RNA_reads_number, spike_in_RNA_s )


timestr = time.strftime("%Y%m%d-%H%M%S")



print(timestr)

newline="\t".join(str(i) for i in("input_file_name","RNA_coverage_normalized", "spike_in_RNA_coverage", "spike_in_RNA_2500391", "time_checked" ))+"\n"
out_put_file.write(newline)
newline="\t".join(str(i) for i in(input_file_name_s,RNA_coverage_normalized, spike_in_RNA_coverage, spike_in_RNA_s, timestr ))+"\n"
out_put_file.write(newline)

out_put_file.close()