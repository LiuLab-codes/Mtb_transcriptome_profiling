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
spike_in_RNA_coverage=1
RNA_coverage=0
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
            if int(position_informaton[0]) <=4411709: RNA_coverage += int(genome[int(position_informaton[0])][3][6])
            if int(position_informaton[0]) > 4411709: spike_in_RNA_coverage +=int(genome[int(position_informaton[0])][3][3])
        
        

        
bed_read_file_1.close()   
spike_in_RNA_coverage=round(spike_in_RNA_coverage/600000,2 )+1
RNA_coverage_normalized= round(RNA_coverage/4456782, 2)
normalization_efficiency= 1


if input_file_name_s=="Rho6_sgRNA__ATC72h_SC3_total_S2": normalization_efficiency=0.823
if input_file_name_s=="Rho8_sgRNA__ATC72h_SC2_total_S3": normalization_efficiency=0.994
if input_file_name_s=="Rho_sg6_ATC_96h_total_S2": normalization_efficiency=1.317
if input_file_name_s=="Rho_sg6_C_96h_total_S1": normalization_efficiency=0.889
if input_file_name_s=="Rho_sg7_ATC_DMSO_total_S6": normalization_efficiency=1.271
if input_file_name_s=="Rho_sg7_ATC_Rifa_total_S7": normalization_efficiency=1.271
if input_file_name_s=="Rho_sg7_ATC_linez_total_S8": normalization_efficiency=0.747
if input_file_name_s=="Rho_sg7_C_DMSO_total_S3": normalization_efficiency=0.954
if input_file_name_s=="Rho_sg7_C_Rifa_total_S4": normalization_efficiency=0.902
if input_file_name_s=="Rho_sg7_C_linez_total_S5": normalization_efficiency=1.194

print(input_file_name_s,RNA_coverage_normalized, normalization_efficiency)

operon_elongation_infor_inputfile= open('operon_elongation_site07172021.txt', 'r')
output_file= open(input_file_name_s+ '_elongation_profiling07282021.txt',  'w')


for line in  operon_elongation_infor_inputfile:
    position_informaton = line.strip().split("\t")
    
    if not position_informaton[1].isdigit():
        position_informaton.pop()
        newline="\t".join(str(i) for i in position_informaton) + "\t"+"\t".join(str(i) for i in("initial_length","initial_expression", "normalized_initial_expression", "normalized_initial_expression_2_spike_RNA",  "elongation_length", "elongation_expression", "normalized_elongation_expression","normalized_elongation_expression_2_spike_RNA", "elongation_ratio",  "operon_elongation_length", "oepron_relative_length", "sense_expression", "sense_2_spike_in_RNA", "antisense_expression","antisense_2_spike_in_RNA", "antisense_ratio", "antisense_ratio_percentage"))+"\n"
        output_file.write(newline)
        
        
        
        
    
    
    
    if not position_informaton[1].isdigit(): continue
    
    
    operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check,leaderless_check, included_gene_number, coding_gene, TSS_peak_site,  tss_extention_end, elongation_check_start_site,	elongation_check_end_site, elongation_check_length, *others =position_informaton
    
    
    #print(operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, included_gene_number, coding_gene, elongation_check_start_site,	elongation_check_end_site, elongation_check_length,)
    
    
    
    if operon_direction=="+":
        
        initial_expression, elongation_expression=(0,0)
        initial_length, elongation_length=(0,0)
        for screen_position in range ( int(TSS_peak_site), int(tss_extention_end)):
            initial_length+=1
            initial_expression+=int(genome[int(screen_position)][3][3])
        
        for screen_position in range ( int(elongation_check_start_site), int(elongation_check_end_site)):
            elongation_length+=1
            elongation_expression+=int(genome[int(screen_position)][3][3])
        
        if elongation_length==0: elongation_length=1
        if initial_length==0: initial_length=1
        initial_expression=round(initial_expression/initial_length,2)
        elongation_expression=round(elongation_expression/elongation_length,2)
        elongation_ratio= round(elongation_expression/(initial_expression+1), 2)
        normalized_initial_expression = round(initial_expression*normalization_efficiency, 2)
        normalized_initial_expression_2_spike_RNA = round(initial_expression/spike_in_RNA_coverage, 2)
        normalized_elongation_expression=round(elongation_expression*normalization_efficiency, 2)
        normalized_elongation_expression_2_spike_RNA=round(elongation_expression/spike_in_RNA_coverage, 2)
        
        highest_coverage=0
        operon_elongation_length=int(operon_end_site)-int(operon_start_site)+1
        
        for position in range(int(operon_start_site), int(elongation_check_start_site)):
            if position<5: continue
            
            expression=int((int(genome[position-2][3][3])+int(genome[position-1][3][3])+int(genome[position-0][3][3])+int(genome[position+1][3][3])+int(genome[position+2][3][3]))/5)
            if expression > highest_coverage: highest_coverage=expression
            
        for position in range(int(elongation_check_start_site), int(operon_end_site)):
            if position<5: continue
            expression=int((int(genome[position-2][3][3])+int(genome[position-1][3][3])+int(genome[position-0][3][3])+int(genome[position+1][3][3])+int(genome[position+2][3][3]))/5)
            
            if expression< highest_coverage*0.33: 
                operon_elongation_length=position-int(operon_start_site)
                
                break
            
        
        
        oepron_relative_length= round(operon_elongation_length/(int(operon_end_site)-int(operon_start_site)+1),2)
        
        
        
        sense_transcription, antisense_transcription=(0,0)
        
        for position in range(int(operon_start_site), int(operon_end_site)):
            sense_transcription+= int(genome[position][3][3])
            antisense_transcription += int(genome[position][3][6])
        sense_2_spike_in_RNA= round((sense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_2_spike_in_RNA= round((antisense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_ratio_percentage=round(antisense_transcription/(sense_transcription+antisense_transcription+1)*100,1)
        sense_transcription= round((sense_transcription/int(operon_length))*normalization_efficiency, 2)
        antisense_transcription = round((antisense_transcription/int(operon_length))*normalization_efficiency, 3)
        
        
        antisense_ratio=0
        if sense_transcription>0: antisense_ratio=round(antisense_transcription/sense_transcription, 3)
        
        
        
            
            
            
        
        
        
        
        #print(operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, included_gene_number, coding_gene, TSS_peak_site,tss_extention_end, elongation_check_start_site,	elongation_check_end_site,initial_length,initial_expression, normalized_initial_expression,  elongation_length, elongation_expression, elongation_ratio,normalized_elongation_expression, operon_elongation_length )
        
        newline="\t".join(str(i) for i in  (operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, leaderless_check, included_gene_number, coding_gene, TSS_peak_site, tss_extention_end, elongation_check_start_site,	elongation_check_end_site,elongation_check_length, initial_length,initial_expression, normalized_initial_expression, normalized_initial_expression_2_spike_RNA,  elongation_length, elongation_expression, normalized_elongation_expression, normalized_elongation_expression_2_spike_RNA,  elongation_ratio, operon_elongation_length, oepron_relative_length, sense_transcription, sense_2_spike_in_RNA, antisense_transcription, antisense_2_spike_in_RNA, antisense_ratio, antisense_ratio_percentage))+"\n"
        
        output_file.write(newline)
        
        
    if operon_direction=="-":
        
        initial_expression, elongation_expression=(0,0)
        initial_length, elongation_length=(0,0)
        for screen_position in range ( int(TSS_peak_site), int(tss_extention_end), -1):
            initial_length+=1
            initial_expression+=int(genome[int(screen_position)][3][6])
        
        for screen_position in range ( int(elongation_check_start_site), int(elongation_check_end_site), -1):
            elongation_length+=1
            elongation_expression+=int(genome[int(screen_position)][3][6])
        
        if elongation_length==0: elongation_length=1
        if initial_length==0: initial_length=1
        initial_expression=round(initial_expression/initial_length,2)
        elongation_expression=round(elongation_expression/elongation_length,2)
        elongation_ratio= round(elongation_expression/(initial_expression+1), 2)
        normalized_initial_expression = round(initial_expression*normalization_efficiency, 2)
        normalized_initial_expression_2_spike_RNA = round(initial_expression/spike_in_RNA_coverage, 2)
        normalized_elongation_expression=round(elongation_expression*normalization_efficiency, 2)
        normalized_elongation_expression_2_spike_RNA=round(elongation_expression/spike_in_RNA_coverage, 2)
        
        
        
        
        highest_coverage=0
        operon_elongation_length=int(operon_end_site)-int(operon_start_site)+1
        
        for position in range(int(operon_end_site), int(elongation_check_start_site), -1):
            
            expression=int((int(genome[position-2][3][6])+int(genome[position-1][3][6])+int(genome[position-0][3][6])+int(genome[position+1][3][6])+int(genome[position+2][3][6]))/5)
            if expression > highest_coverage: highest_coverage=expression
            
        for position in range(int(elongation_check_start_site), int(operon_start_site), -1):
            
            expression=int((int(genome[position-2][3][6])+int(genome[position-1][3][6])+int(genome[position-0][3][6])+int(genome[position+1][3][6])+int(genome[position+2][3][6]))/5)
            
            if expression< highest_coverage*0.33: 
                operon_elongation_length=int(operon_end_site)-position
                
                break
        
        
        oepron_relative_length= round(operon_elongation_length/(int(operon_end_site)-int(operon_start_site)+1),2)
        
        
        sense_transcription, antisense_transcription=(0,0)
        
        for position in range(int(operon_start_site), int(operon_end_site)):
            sense_transcription+= int(genome[position][3][6])
            antisense_transcription += int(genome[position][3][3])
        sense_2_spike_in_RNA= round((sense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_2_spike_in_RNA= round((antisense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_ratio_percentage=round(antisense_transcription/(sense_transcription+antisense_transcription+1)*100,1) 
        sense_transcription= round((sense_transcription/int(operon_length))*normalization_efficiency, 3) 
        #sense_transcription_2_spike_RNA= round((sense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_transcription = round((antisense_transcription/int(operon_length))*normalization_efficiency, 3)
        #antisense_transcription_2_spike_RNA = round((antisense_transcription/int(operon_length))/spike_in_RNA_coverage, 3)
        antisense_ratio=0
        if sense_transcription>0: antisense_ratio=round(antisense_transcription/sense_transcription, 3)
        
        
        #print(operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, included_gene_number, coding_gene, TSS_peak_site,tss_extention_end, elongation_check_start_site,	elongation_check_end_site,initial_length,initial_expression, normalized_initial_expression,  elongation_length, elongation_expression, elongation_ratio,normalized_elongation_expression )
        
        newline="\t".join(str(i) for i in  (operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, leaderless_check, included_gene_number, coding_gene, TSS_peak_site,tss_extention_end,  elongation_check_start_site,	elongation_check_end_site,elongation_check_length, initial_length,initial_expression, normalized_initial_expression,  normalized_initial_expression_2_spike_RNA, elongation_length, elongation_expression, normalized_elongation_expression, normalized_elongation_expression_2_spike_RNA, elongation_ratio, operon_elongation_length, oepron_relative_length, sense_transcription, sense_2_spike_in_RNA, antisense_transcription,antisense_2_spike_in_RNA,  antisense_ratio, antisense_ratio_percentage))+"\n"
        
        output_file.write(newline)
    
        
        
            
        
        
        












