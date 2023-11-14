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
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

TSS_transcript_infor = {}
genome = {}
genome_seq_signal={}
#cwd = os.getcwd()
#input_file_name = os.path.splitext(str(sys.argv[1]))[0]
#pattern = input_file_name.split("_PE")
#input_file_name=pattern[0]

#print(input_file_name)







#home_directory = os.getenv("HOME")
for num in range (0,4471711):
#for num in range (0,471711):
    #genome[num] = [0,0,0,0,0,0,0,0]
    genome_seq_signal[num] = [num,[], [], []]
    if num > 210000: break
    
print ('finish the first step of data initialization\n')



#output_file=open(input_file_name+"_operon_coverage.txt", 'w')
#bedfile_full_directory = cwd+'/'+str(sys.argv[1])
#print ("start to read the bedpe file: " + bedfile_full_directory + "\n")        
bed_input_file = open("S05_total_RD_TB_WT_log_bed_reading.txt", 'r')
#1st file
group1_name ="WT_total"
group1_name_total_coverage =0
for line_bed in  bed_input_file:
    position_informaton_bed = line_bed.strip().split("\t")
        
    if (len(position_informaton_bed) >4) :
        
        position, *coverage_infor = position_informaton_bed
        if int(position)>210000: break
        if int(position) > 4471710 : continue  #pay attention to this restriction
        genome_seq_signal[int(position)][1]=coverage_infor
        #print(genome_seq_signal[int(position)][1])
        #group1_name_total_coverage += (int(coverage_infor[2])+int(coverage_infor[5]))
bed_input_file.close
# format of input_bed_data: 67	[0	0	25	0	0	19]



number=1
operon_infor_input_file = open('operon_annoation07032021.txt', "r")
for line_infor in  operon_infor_input_file:
    position_informaton = line_infor.strip().split("\t")
    if len(position_informaton)<6 : continue
    if not position_informaton[1].isdigit(): continue
    
    operon_name, original_operon_start_site, original_operon_end_site, operon_length, operon_direction=(position_informaton[0],int(position_informaton[1]),int(position_informaton[2]),int(position_informaton[3]),position_informaton[4] )
    
    
    
    if (original_operon_end_site-original_operon_start_site+1)<300: continue
    if original_operon_end_site > 210000: break
    if operon_direction == "+":
        
        max_coverage=1
        
        for  screen_position in range (original_operon_start_site, (original_operon_end_site+1)):
            
            if int(genome_seq_signal[int(screen_position)][1][2])>max_coverage: max_coverage=int(genome_seq_signal[int(screen_position)][1][2])
            
            #if screen_position>original_operon_start_site+600 or screen_position>original_operon_start_site+operon_length*0.6:
                #break
        initial_expression=1        
        for   screen_position in range (original_operon_start_site, (original_operon_start_site+200)):
            initial_expression+=int(genome_seq_signal[int(screen_position)][1][2])
        initial_expression=int(initial_expression/200)
            
        if initial_expression<50: continue
        
        print(operon_name, original_operon_start_site, original_operon_end_site, operon_length, max_coverage)
    
        X_nucletide_posiiton = []
        Y_1_coverage_p = []
        Y_1_coverage_n = []
        G1_coverage_origin_p=[]
        G1_coverage_origin_n=[]
        
        gene_cover_start_left=original_operon_start_site
        gene_cover_start_right=original_operon_end_site
        
        #print(gene_cover_start_left,gene_cover_start_right)
        for screen_position in range (gene_cover_start_left, (gene_cover_start_right-10)):
                G1_coverage_origin_p.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][1][2])/max_coverage*100])
                #G1_coverage_origin_n.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][1][6])/max_coverage*100])
                #G2_coverage_origin_p.append([(int(screen_position)), int(genome[int(screen_position)][4][3])])
                #G2_coverage_origin_n.append([(int(screen_position)), int(genome[int(screen_position)][4][6])])
                #G3_coverage_origin_p.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][2])/group3_name_total_coverage*1000000000])
                #G3_coverage_origin_n.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][5])/group3_name_total_coverage*1000000000])
        
        
        
        
        G1_coverage_p_array = numpy.array(G1_coverage_origin_p)
        #G1_coverage_n_array = numpy.array(G1_coverage_origin_n)
        
        
        
        
        
        
        
        
        
        
        
        
        x_Y_1_coverage_p,y_Y_1_coverage_p=numpy.array(G1_coverage_p_array).T
        
        
        
        plt.plot(x_Y_1_coverage_p, y_Y_1_coverage_p,  color="pink", linewidth=0.1)
        
    
    
    
    
plt.ylim((0,120))   
plt.xlim((0,3000))
plt.ylabel('Normalized_coverage')
    #plt.xlabel( 'genome position of ' + operon_name +"(" + operon_direction + ")")
plt.savefig("total_p_coverage_operon"+".png", dpi=300)
    
        
        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

