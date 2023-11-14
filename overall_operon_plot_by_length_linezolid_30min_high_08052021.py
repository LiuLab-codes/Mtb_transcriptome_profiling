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
    genome_seq_signal[num] = [num,"", "", ""]
    #if num > 210000: break
    
print ('finish the first step of data initialization\n')



#output_file=open(input_file_name+"_operon_coverage.txt", 'w')
#bedfile_full_directory = cwd+'/'+str(sys.argv[1])
#print ("start to read the bedpe file: " + bedfile_full_directory + "\n")        





bed_read_file_1 = open('BSL2TB_DMSO_10min_total_S9_bed_reading.txt', 'r')

for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        position, *coverage_infor = position_informaton
        #if int(position)>210000: break
        if int(position) > 4471710 : continue
        #print(position_informaton)
        #if int(position_informaton[0]) > 210000: break
        genome_seq_signal[int(position_informaton[0])][1] = coverage_infor

        
bed_read_file_1.close()   


bed_read_file_2 = open('linezolid_30min_total_S4_bed_reading.txt', 'r')

for line in  bed_read_file_2:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        position, *coverage_infor = position_informaton
        #if int(position)>210000: break
        if int(position) > 4471710 : continue
        #if int(position_informaton[0]) > 210000: break
        genome_seq_signal[int(position_informaton[0])][2] = coverage_infor
        
bed_read_file_2.close() 







average_3000={}
average_3000_2nd={}
for site in range(0, 3000):
    average_3000[site]=[[],[]]
    average_3000_2nd[site]=[[],[]]

efficiency_file1= 0.797
efficiency_file2= 1.02

number=1
operon_infor_input_file = open('operon_annoation07112021.txt', "r")
for line_infor in  operon_infor_input_file:
    position_informaton = line_infor.strip().split("\t")
    if len(position_informaton)<6 : continue
    if not position_informaton[1].isdigit(): continue
    
    operon_name, original_operon_start_site, original_operon_end_site, operon_length, operon_direction=(position_informaton[0],int(position_informaton[1]),int(position_informaton[2]),int(position_informaton[3]),position_informaton[4] )
    
    
    
    
    
    initial_expression, elongation_ratio= (float(position_informaton[10]), float(position_informaton[11]))
    
    
    if (original_operon_end_site-original_operon_start_site+1)<500: continue
    #if original_operon_end_site > 210000: break
    
    if initial_expression<20: continue
    
    if elongation_ratio <=0.3: continue
    
  
    if operon_direction == "+":
        
        max_coverage=1
        max_coverage_2nd=1
        
        for  screen_position in range (original_operon_start_site, (original_operon_end_site+1)):
            
            if int(genome_seq_signal[int(screen_position)][1][2])*efficiency_file1>max_coverage: max_coverage=int(genome_seq_signal[int(screen_position)][1][2])*efficiency_file1
            if int(genome_seq_signal[int(screen_position)][2][2])*efficiency_file2>max_coverage: max_coverage=int(genome_seq_signal[int(screen_position)][2][2])*efficiency_file2
            
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
        #G1_coverage_origin_n=[]
        G2_coverage_origin_p=[]
        #G2_coverage_origin_n=[]
        gene_cover_start_left=original_operon_start_site
        gene_cover_start_right=original_operon_end_site
        
        #print(gene_cover_start_left,gene_cover_start_right)
        for screen_position in range (gene_cover_start_left, (gene_cover_start_right-10)):
                G1_coverage_origin_p.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][1][2])/max_coverage*100*efficiency_file1])
                G2_coverage_origin_p.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][2][2])/max_coverage*100*efficiency_file2])
                #G1_coverage_origin_n.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][1][6])/max_coverage*100])
                #G2_coverage_origin_p.append([(int(screen_position)), int(genome[int(screen_position)][4][3])])
                #G2_coverage_origin_n.append([(int(screen_position)), int(genome[int(screen_position)][4][6])])
                #G3_coverage_origin_p.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][2])/group3_name_total_coverage*1000000000])
                #G3_coverage_origin_n.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][5])/group3_name_total_coverage*1000000000])
                if operon_length >= 500:
                     
                    if (int(screen_position)-gene_cover_start_left) <3000:
                        average_3000[(int(screen_position)-gene_cover_start_left)][1].append(int(genome_seq_signal[int(screen_position)][1][2])*efficiency_file1/max_coverage*100)
                        average_3000_2nd[(int(screen_position)-gene_cover_start_left)][1].append(int(genome_seq_signal[int(screen_position)][2][2])*efficiency_file2/max_coverage*100)
                        
                    
        
        
        
        
        G1_coverage_p_array = numpy.array(G1_coverage_origin_p)
        G2_coverage_p_array = numpy.array(G2_coverage_origin_p)
        #G1_coverage_n_array = numpy.array(G1_coverage_origin_n)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        x_Y_1_coverage_p,y_Y_1_coverage_p=numpy.array(G1_coverage_p_array).T
        x_Y_2_coverage_p,y_Y_2_coverage_p=numpy.array(G2_coverage_p_array).T
        
        if operon_length >= 500: 
            
            line_color="red"
            line_color_2nd="blue"
            
            
        
        #plt.plot(x_Y_1_coverage_p, y_Y_1_coverage_p,  color=line_color, linewidth=0.02, alpha=0.5)
        #plt.plot(x_Y_2_coverage_p, y_Y_2_coverage_p,  color=line_color_2nd, linewidth=0.02, alpha=0.5)
        number+=1
#print(average_1500)       

    if operon_direction == "-":
        
        max_coverage=1
        max_coverage_2nd=1
        for  screen_position in range (original_operon_end_site, (original_operon_start_site+1), -1):
            
            if int(genome_seq_signal[int(screen_position)][1][5])*efficiency_file1>max_coverage: max_coverage=int(genome_seq_signal[int(screen_position)][1][5])*efficiency_file1
            if int(genome_seq_signal[int(screen_position)][2][5])*efficiency_file2>max_coverage_2nd: max_coverage_2nd=int(genome_seq_signal[int(screen_position)][2][5])*efficiency_file2
            
            #if screen_position<original_operon_end_site-600 or screen_position<original_operon_end_site-operon_length*0.6:
                #break
        initial_expression=1        
        for   screen_position in range (original_operon_end_site, (original_operon_end_site-200), -1):
            initial_expression+=int(genome_seq_signal[int(screen_position)][1][5])
        initial_expression=int(initial_expression/200)
            
        if initial_expression<50: continue
        
        print(operon_name, original_operon_start_site, original_operon_end_site, operon_length, max_coverage)
    
        X_nucletide_posiiton = []
        Y_1_coverage_p = []
        Y_1_coverage_n = []
        #G1_coverage_origin_p=[]
        G1_coverage_origin_n=[]
        G2_coverage_origin_n=[]
        gene_cover_start_left=original_operon_start_site
        gene_cover_start_right=original_operon_end_site
        
        #print(gene_cover_start_left,gene_cover_start_right)
        for screen_position in range (gene_cover_start_right, (gene_cover_start_left+10), -1):
                G1_coverage_origin_n.append([(gene_cover_start_right-int(screen_position)), int(genome_seq_signal[int(screen_position)][1][5])/max_coverage*100*efficiency_file1])
                G2_coverage_origin_n.append([(gene_cover_start_right-int(screen_position)), int(genome_seq_signal[int(screen_position)][2][5])/max_coverage_2nd*100*efficiency_file2])
                #G1_coverage_origin_n.append([(int(screen_position)-gene_cover_start_left), int(genome_seq_signal[int(screen_position)][1][6])/max_coverage*100])
                #G2_coverage_origin_p.append([(int(screen_position)), int(genome[int(screen_position)][4][3])])
                #G2_coverage_origin_n.append([(int(screen_position)), int(genome[int(screen_position)][4][6])])
                #G3_coverage_origin_p.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][2])/group3_name_total_coverage*1000000000])
                #G3_coverage_origin_n.append([(int(screen_position)), int(genome_seq_signal[int(screen_position)][3][0][5])/group3_name_total_coverage*1000000000])
                if operon_length >= 500: 
                    
                    if (gene_cover_start_right-int(screen_position)) <3000:
                        average_3000[(gene_cover_start_right-int(screen_position))][1].append(int(genome_seq_signal[int(screen_position)][1][5])/max_coverage*100*efficiency_file1)
                        average_3000_2nd[(gene_cover_start_right-int(screen_position))][1].append(int(genome_seq_signal[int(screen_position)][2][5])/max_coverage_2nd*100*efficiency_file2)
                    
        
        
        
        
        #G1_coverage_p_array = numpy.array(G1_coverage_origin_p)
        G1_coverage_n_array = numpy.array(G1_coverage_origin_n)
        G2_coverage_n_array = numpy.array(G2_coverage_origin_n)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        x_Y_1_coverage_n,y_Y_1_coverage_n=numpy.array(G1_coverage_n_array).T
        x_Y_2_coverage_n,y_Y_2_coverage_n=numpy.array(G2_coverage_n_array).T
        
        
        if operon_length > 500: 
            
            line_color="orange"
            line_color_2nd="blue"
        
        #plt.plot(x_Y_1_coverage_n, y_Y_1_coverage_n,  color=line_color, linewidth=0.02, alpha=0.5)
        #plt.plot(x_Y_2_coverage_n, y_Y_2_coverage_n,  color=line_color, linewidth=0.02, alpha=0.5)
        number+=1
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
all_average_3000=[]
for site in range(0, 3000):
    average_coverage_3000=0
    if len(average_3000[site][1])>0: average_coverage_3000 =   numpy.median(average_3000[site][1])        
    
    error_3000=0
    if len(average_3000[site][1])>2:error_3000=statistics.stdev(average_3000[site][1])
    
    all_average_3000.append([site, average_coverage_3000, error_3000 ])
    

all_average_3000_array=numpy.array(all_average_3000)

x_3000, y_3000, error_3000 = numpy.array(all_average_3000_array).T

plt.plot(x_3000, y_3000,  color="red", label="DMSO", linewidth=1.5)
plt.fill_between(x_3000, y_3000-error_3000,y_3000+error_3000, color="coral", alpha=0.25 )
#plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)    




all_average_3000_2nd=[]
for site in range(0, 3000):
    average_coverage_3000_2nd=0
    if len(average_3000_2nd[site][1])>0: average_coverage_3000_2nd =   numpy.median(average_3000_2nd[site][1])         
    
    error_3000_2nd=0
    if len(average_3000_2nd[site][1])>2:error_3000_2nd=statistics.stdev(average_3000_2nd[site][1])
    
    all_average_3000_2nd.append([site, average_coverage_3000_2nd, error_3000_2nd ])
    

all_average_3000_2nd_array=numpy.array(all_average_3000_2nd)

x_3000_2nd, y_3000_2nd, error_3000_2nd = numpy.array(all_average_3000_2nd_array).T

plt.plot(x_3000_2nd, y_3000_2nd,  color="blue", label="Linezolid", linewidth=1.5)
plt.fill_between(x_3000_2nd, y_3000_2nd-error_3000_2nd,y_3000_2nd+error_3000_2nd, color="lightblue", alpha=0.25  )
    
plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)  





    
plt.ylim((0,100))   
plt.xlim((0,3000))
plt.ylabel('Normalized_coverage')
plt.xlabel( "Distance to the operon start site (bp)")
    #plt.xlabel( 'genome position of ' + operon_name +"(" + operon_direction + ")")
plt.savefig("overall_linezolid_30min_high"+"_N_"+str(number)+"_all_operon_plot_by_length08062021"+".pdf", dpi=600)
    
        
print("total counted transcritps is ", number )        
    
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

