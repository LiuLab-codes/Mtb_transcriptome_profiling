#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import glob, os

## This script is designed to gather all the TTS identified from samples located in the same directory. The criteria for filtering may be adjusted as needed.




genome = {}


#input_file_name = os.path.splitext(str(sys.argv[1]))[0]
home_directory = os.getenv("HOME")
for num in range (0,4471712):
    genome[num] = ["",0,0,[],[],[],[],[],[]]
    #genome_seq_signal[num] = [num,0,0,0,0,0,0]
    #if num > 220000 : break

print ('finish the first step of data initialization\n')


input_direct_name = os.path.splitext(str(sys.argv[1]))[0]


TTS= {}
os.chdir(input_direct_name)
sample_serial =0

for file in glob.glob("*tts_extraction06292023.txt"):
    sample_serial +=1
    file_name = os.path.splitext(file)[0]
    pattern = file_name.split("_tts")
    input_file_name_s=pattern[0]
    print(input_file_name_s)
    
    input_file= open (file, 'r')
    for line in  input_file:
        informaton = line.strip().split("\t")
        if not informaton[0].isdigit(): continue
        tts_position, tts_direction, tts_intensity, tts_efficiency = (informaton[0],informaton[1],informaton[2], informaton[3])
        if int(tts_intensity)<10: continue
        #print(tts_position, tts_direction, tts_intensity)
        if tts_direction == "+":
            genome[int(tts_position)][1]+=1
            genome[int(tts_position)][3].append(int(tts_intensity))
            genome[int(tts_position)][5].append(float(tts_efficiency))
            genome[int(tts_position)][7].append(input_file_name_s)
        if tts_direction == "-":
            genome[int(tts_position)][2]+=1
            genome[int(tts_position)][4].append(int(tts_intensity))
            genome[int(tts_position)][6].append(float(tts_efficiency))
            genome[int(tts_position)][8].append(input_file_name_s)
out_put= open("TTS_collection06292023.txt",'w')  

newline= "\t".join(str(i) for i in ( "genome_position", "direction", "TTS_intensity","Termination_efficiency","repeat_times",  "source"))
out_put.write(newline + "\n")           
next_check =0
for num in range (0,4471712):
    if num < next_check : continue
    #if genome[num][1] ==0 and genome[num][2]==0: continue   
    if  genome[num][1] !=0:
        peak= num
        average_intensity=0
        total_number =0
        efficiency=0
        source=[]
        for position_shift in range (-5,6):
            
            
            if genome[num+position_shift][1]>0:
                average_intensity+=sum(genome[num+position_shift][3])
                source+=genome[num+position_shift][7]
                total_number += genome[num+position_shift][1]
                efficiency+=sum(genome[num+position_shift][5])
            if sum(genome[num+position_shift][3])> sum(genome[peak][3]):
                peak= num+position_shift
                
                
        average_intensity = int(average_intensity/total_number)
        efficiency=round(efficiency/total_number, 2)
        source_infor= ":".join(str(i) for i in source)
        #peak+=1
        #print(peak, "+", average_intensity, total_number, source_infor)
        next_check = num+5
        newline =  "\t".join(str(i) for i in (peak, "+", average_intensity, efficiency, total_number, source_infor)) +"\n"
        if total_number >= 2: out_put.write(newline)
    
next_check =0
for num in range (0,4471712):
    if num < next_check : continue
    #if genome[num][1] ==0 and genome[num][2]==0: continue   
    if  genome[num][2] !=0:
        peak= num
        average_intensity=0
        total_number =0
        efficiency=0
        source=[]
        for position_shift in range (-5,6):
            
            
            if genome[num+position_shift][2]>0:
                average_intensity+=sum(genome[num+position_shift][4])
                source+=genome[num+position_shift][8]
                efficiency+=sum(genome[num+position_shift][6])
                total_number += genome[num+position_shift][2]
            if sum(genome[num+position_shift][4])> sum(genome[peak][4]):
                peak= num+position_shift
                
                
        average_intensity = int(average_intensity/total_number)
        efficiency=round(efficiency/total_number, 2)
        source_infor= ":".join(str(i) for i in source)    
        
        
        #print(peak, "-", average_intensity, total_number, source_infor)
        next_check = num+5
        
        newline =  "\t".join(str(i) for i in (peak, "-", average_intensity, efficiency,total_number, source_infor))+"\n"
        if total_number >= 2: out_put.write(newline)
out_put.close()