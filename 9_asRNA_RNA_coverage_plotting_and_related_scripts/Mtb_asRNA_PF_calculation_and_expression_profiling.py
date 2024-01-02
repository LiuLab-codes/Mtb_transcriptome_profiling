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

home_directory = os.getenv("HOME")
for num in range (0,4471711):

    genome[num] = ["",0,0,"","","",[],[]]
    
    
print ('finish the first step of data initialization\n')




TSS_input_file =  open('Mtb_TSS_final.txt', 'r')

for line in  TSS_input_file:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: continue
        #if int(position_informaton[0]) >21000:
            #continue
    
        if (len(position_informaton) >1):
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
RNA_coverage=0
bed_read_file_1 = open(input_file_name, 'r')
for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        genome[int(position_informaton[0])][3] = position_informaton
        
        if   1471585 < int(position_informaton[0]) < 1477050 : continue
        
        RNA_coverage += int(genome[int(position_informaton[0])][3][3])
        RNA_coverage += int(genome[int(position_informaton[0])][3][6])
        
        
        

        
bed_read_file_1.close()   

RNA_coverage_normalized= round(RNA_coverage/4462246, 2)

print(input_file_name_s,RNA_coverage_normalized)


as_infor_inputfile= open('Mtb_asRNA_information_log_phase.txt', 'r')
output_file= open(input_file_name_s+ '_as_expression_profiling_PF_calculation.txt',  'w')
newline="\t".join(str(i) for i in  ( "TU_start_site", "TU_end_site","TU_direction", "AS_transcript_start", "AS_transcript_end", "AS_transcript_direction",  "AS_3end_information","AS_length",  "as_max_coverage", "as_coverage", "actual_length","start_site_expression", "as_tss_p", "upstream_coverage_average", "down_coverage_average",  "as_PF" ))+"\n"

output_file.write(newline)

for as_line in as_infor_inputfile:
    
    position_informaton = as_line.strip().split("\t")
    
    if not position_informaton[1].isdigit():continue
    TU_start_site, TU_end_site, TU_direction, AS_transcript_start, AS_transcript_end, AS_transcript_direction, AS_3end_information= (position_informaton[0],position_informaton[1],position_informaton[2],int(position_informaton[3]),int(position_informaton[4]),position_informaton[5], position_informaton[6])
    
    if AS_transcript_direction =="+":
        
        as_coverage, as_length= (0, int(position_informaton[7]))
        
        actual_length, as_max_coverage, start_site_expression, check_length=(as_length-100,0,0,0)
        for screen_position in range ( int(AS_transcript_start), int(AS_transcript_end)):
            if screen_position<1 or screen_position> 4471710 :continue
            as_coverage+=int(genome[int(screen_position)][3][3])
            if as_max_coverage < int(genome[int(screen_position)][3][3]):as_max_coverage= int(genome[int(screen_position)][3][3])
            check_length+=1
        as_coverage=int(as_coverage/check_length)
        normalzied_as_coverage=round(as_coverage/RNA_coverage_normalized, 2)
        for screen_position in range ( int(AS_transcript_start), int(AS_transcript_start)+100):
            if screen_position<1 or screen_position> 4471710 :continue
            start_site_expression+=int(genome[int(screen_position)][3][3])
        start_site_expression=int(start_site_expression/100)
            
        
        
        if start_site_expression>20:
            
            for screen_position in range ( int(AS_transcript_end)-100, int(AS_transcript_end)+100):
                if screen_position<1 or screen_position> 4471710 :continue
                actual_length+=1
                
                if int(genome[int(screen_position)][3][3]) < 0.25*as_max_coverage:
                    break
                    
        as_middle_coverage_site=int(0.5*AS_transcript_start+0.5*AS_transcript_end)
        if start_site_expression>20:
            for screen_position in range ( int(AS_transcript_start)+50, int(AS_transcript_start)+actual_length):
                if screen_position<1 or screen_position> 4471710 :continue
                as_middle_coverage_site=screen_position
                
                if int(genome[int(screen_position)][3][3]) < 0.5*start_site_expression:
                    break
            
                    
                  
        
        
        (upstream_length, upstream_coverage, downstream_length,dowstream_coverage)=(0,1,1,1)
        
        for screen_position in range ( AS_transcript_start,AS_transcript_start+200):
            upstream_length+=1
            upstream_coverage+=int(genome[int(screen_position)][3][3])
        upstream_coverage_average= upstream_coverage/upstream_length
        
        
        
        as_tss_p=0
        for screen_position in range (AS_transcript_start+10, AS_transcript_start+800):
            
            if genome[screen_position][1]>0: as_tss_p+=1
        
        if as_tss_p==0:
            
            for screen_position in range ( AS_transcript_start+500,AS_transcript_start+700):
                downstream_length+=1
                dowstream_coverage+=int(genome[int(screen_position)][3][3])
            
            
        
        
       
        down_coverage_average= dowstream_coverage/downstream_length
        
        
        as_PF= round(down_coverage_average/upstream_coverage_average,2)
        
        
        
        
        
        print( TU_start_site, TU_end_site, TU_direction, AS_transcript_start, AS_transcript_end, AS_transcript_direction, AS_3end_information, as_length, as_max_coverage, as_coverage,normalzied_as_coverage, actual_length, start_site_expression, as_tss_p, upstream_coverage_average, down_coverage_average,  as_PF)
        
        newline="\t".join(str(i) for i in (TU_start_site, TU_end_site, TU_direction, AS_transcript_start, AS_transcript_end, AS_transcript_direction, AS_3end_information, as_length, as_max_coverage, as_coverage, actual_length, start_site_expression, as_tss_p, upstream_coverage_average, down_coverage_average,  as_PF))+"\t"+input_file_name_s+"\n"
    
        output_file.write(newline)
        
        
        
    if AS_transcript_direction =="-":
        
        as_coverage_n, as_length_n= (0, int(position_informaton[7]))
        
        actual_length_n, as_max_coverage_n, start_site_expression_n, check_length_n=(as_length_n-100,0,0,0)
        for screen_position in range ( int(AS_transcript_end), int(AS_transcript_start), -1):
            if screen_position<1 or screen_position> 4471710 :continue
            as_coverage_n+=int(genome[int(screen_position)][3][6])
            if as_max_coverage_n < int(genome[int(screen_position)][3][6]):as_max_coverage_n= int(genome[int(screen_position)][3][6])
            check_length_n+=1
        as_coverage_n=int(as_coverage_n/check_length_n)
        normalzied_as_coverage_n=round(as_coverage_n/RNA_coverage_normalized, 2)
        
        for screen_position in range ( int(AS_transcript_end),  int(AS_transcript_end)-100, -1):
            if screen_position<1 or screen_position> 4471710 :continue
            start_site_expression_n+=int(genome[int(screen_position)][3][6])
        start_site_expression_n=int(start_site_expression_n/100)
            
        
        
        if start_site_expression_n>20:
            
            for screen_position in range ( int(AS_transcript_start)+100, int(AS_transcript_start)-100, -1):
                if screen_position<1 or screen_position> 4471710 :continue
                actual_length_n+=1
                
                if int(genome[int(screen_position)][3][6]) < 0.25*as_max_coverage_n:
                    break
                    
        
        
       
        
        
        (upstream_length, upstream_coverage, downstream_length,dowstream_coverage)=(0,1,1,1)
        
        for screen_position in range ( AS_transcript_end,AS_transcript_end-200, -1):
            upstream_length+=1
            upstream_coverage+=int(genome[int(screen_position)][3][6])
        upstream_coverage_average= upstream_coverage/upstream_length
        
        
        
        as_tss_n=0
        for screen_position in range (AS_transcript_end-10, AS_transcript_end-800, -1):
            if screen_position<1 or screen_position> 4471710 :continue
            if genome[screen_position][2]>0: as_tss_n+=1
        
        if as_tss_n==0:
            
            for screen_position in range ( AS_transcript_end-500,AS_transcript_end-700, -1):
                downstream_length+=1
                dowstream_coverage+=int(genome[int(screen_position)][3][6])
            
            
        
        
       
        down_coverage_average= dowstream_coverage/downstream_length
        
        
        as_PF= round(down_coverage_average/upstream_coverage_average,2)
        
        
        
        
        
        print( TU_start_site, TU_end_site, TU_direction, AS_transcript_start, AS_transcript_end, AS_transcript_direction,  AS_3end_information, as_length_n, as_max_coverage_n, as_coverage_n,actual_length_n,  start_site_expression_n, as_tss_n, upstream_coverage_average,down_coverage_average, as_PF)
        
        newline="\t".join(str(i) for i in ( TU_start_site, TU_end_site, TU_direction, AS_transcript_start, AS_transcript_end, AS_transcript_direction, AS_3end_information, as_length_n, as_max_coverage_n, as_coverage_n,actual_length_n, start_site_expression_n, as_tss_n, upstream_coverage_average,down_coverage_average, as_PF) )+"\t"+input_file_name_s+"\n"
    
        output_file.write(newline)
                
                
                
                
        
            
            
            
        
        
        
        
    
    
    
    
    
    


















