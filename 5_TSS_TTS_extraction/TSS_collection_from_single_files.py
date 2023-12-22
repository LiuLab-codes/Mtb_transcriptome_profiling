#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import glob, os



# This script is designed to gather all the TSS identified from samples located in the same directory. The criteria for filtering may be adjusted as needed.

genome = {}



home_directory = os.getenv("HOME")
for num in range (0,4471712):
    genome[num] = ["",0,0,[],[],[],[]]
    #genome_seq_signal[num] = [num,0,0,0,0,0,0]
    #if num > 220000 : break

print ('finish the first step of data initialization\n')


input_direct_name = os.path.splitext(str(sys.argv[1]))[0]


TSS= {}
os.chdir(input_direct_name)
sample_serial =0

for file in glob.glob("*tss_extraction07282023.txt"):
    sample_serial +=1
    file_name = os.path.splitext(file)[0]
    pattern = file_name.split("_tss_extraction")
    input_file_name_s=pattern[0]
    print(input_file_name_s)
    
    input_file= open (file, 'r')
    for line in  input_file:
        informaton = line.strip().split("\t")
        if informaton[2]=="TSS_intensity": continue
        tss_position, tss_direction, tss_intensity = (informaton[0],informaton[1],informaton[2])
        if int(tss_intensity)<10: continue
        #print(tss_position, tss_direction, tss_intensity)
        if tss_direction == "+":
            genome[int(tss_position)][1]+=1
            genome[int(tss_position)][3].append(int(tss_intensity))
            genome[int(tss_position)][5].append(input_file_name_s)
            
        if tss_direction == "-":
            genome[int(tss_position)][2]+=1
            genome[int(tss_position)][4].append(int(tss_intensity))
            genome[int(tss_position)][6].append(input_file_name_s)
            
out_put= open("TSS_collection07282023.txt",'w')             
next_check =0
previous_TSS=-100
for num in range (0,4411406):
    if num < next_check : continue
    #if genome[num][1] ==0 and genome[num][2]==0: continue   
    if  genome[num][1] !=0:
        peak= num
        average_intensity=0
        total_number =0
        source=[]
        for position_shift in range (0,13):
            
            
            if genome[num+position_shift][1]>0:
                average_intensity+=sum(genome[num+position_shift][3])
                source+=genome[num+position_shift][5]
                total_number += genome[num+position_shift][1]
            if sum(genome[num+position_shift][3])> sum(genome[peak][3]):
                peak= num+position_shift
                
                
        average_intensity = int(average_intensity/total_number)
            
        source_infor= ":".join(str(i) for i in source)
        distance_to_previous_positive_tss= peak-previous_TSS
        #peak+=1
        #print(peak, "+", average_intensity, total_number, source_infor)
        next_check = num+13
        newline =  "\t".join(str(i) for i in (peak, "+", average_intensity, total_number, distance_to_previous_positive_tss, source_infor)) +"\n"
        if total_number >= 2:       # The number of replicates may vary as required.
            out_put.write(newline)
            previous_TSS=peak
        
    
next_check =4411406
previous_TSS_n= 4411406    
for num in range (4411406,0,-1):
    if num > next_check : continue
    #if genome[num][1] ==0 and genome[num][2]==0: continue   
    if  genome[num][2] !=0:
        peak= num
        average_intensity=0
        total_number =0
        source=[]
        for position_shift in range (0,13):
            
            
            if genome[num-position_shift][2]>0:
                average_intensity+=sum(genome[num-position_shift][4])
                source+=genome[num-position_shift][6]
                total_number += genome[num-position_shift][2]
            if sum(genome[num-position_shift][4])> sum(genome[peak][4]):
                peak= num-position_shift
                
                
        average_intensity = int(average_intensity/total_number)
        source_infor= ":".join(str(i) for i in source)    
        distance_to_previous_negative_tss= previous_TSS_n-peak
        
        #print(peak, "-", average_intensity, total_number, source_infor)
        next_check = num-13
        
        newline =  "\t".join(str(i) for i in (peak, "-", average_intensity, total_number, distance_to_previous_negative_tss, source_infor))+"\n"
        if total_number >= 2: # The number of replicates may vary as required.
            
            out_put.write(newline)
            previous_TSS_n=peak
            
out_put.close()