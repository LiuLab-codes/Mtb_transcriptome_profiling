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
    if int(position_informaton[0]) > 210000 : break
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


TSS_input_file =  open('TSS_for_operon08052021.txt', 'r')

for line in  TSS_input_file:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: continue
        #if int(position_informaton[0]) >21000:
            #continue
    
        if (len(position_informaton) >=3):
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
        
bed_read_file_1.close()        

RNA_abundance_around_selected_tss={}
for site in range(0,1001):
    RNA_abundance_around_selected_tss[site]=[[],[]]
    
number=1

tss_select_file  = open("select_tss_for_RNAP_analysis04112021.txt", 'r')

for tss_line in tss_select_file:
    tss_information=tss_line.strip().split("\t")
    tss_start_site, tss_strand, tss_intensity, tss_category, upstream_positive_tss_distance, upstream_negative_tss_distance, downstream_positive_tss_distance, downstream_negative_tss_distance=tss_information
    
    #if tss_category=="gene_antisense": continue
    
    tss_start_site=int(tss_start_site)

    
    if tss_strand == "+":
        
        abudance_around_tss=0
    
        for position_shift in range(0, 100):
        
            if tss_start_site+position_shift <=0 or tss_start_site+position_shift > 4471710: continue
            abudance_around_tss+=int(genome[(tss_start_site+position_shift)][3][3])
        
        
        if  abudance_around_tss <3000: continue
        if  abudance_around_tss >5000000: continue
        
        abudance_around_tss=abudance_around_tss/100
        
        check_end = tss_start_site+1000
        
        
        for position_shift in range((tss_start_site+100), tss_start_site+1000):
    
            if genome[position_shift][1]>50:
        
                upstream_tss_coverage, downstream_tss_coverage=(0,1)
        
                for position_shift_tss_check in range (0, 20):
            
                    upstream_tss_coverage +=int(genome[position_shift-position_shift_tss_check][3][3])
                    downstream_tss_coverage +=int(genome[position_shift+position_shift_tss_check+1][3][3])
                if downstream_tss_coverage> upstream_tss_coverage*2 :
                    if position_shift < check_end:
                
                        check_end=position_shift-10
                        break
                        
        for position_shift in range((tss_start_site), (check_end+1)):
            RNA_coverage_site= int(genome[position_shift][3][3])
            #if RNA_coverage_site >5000: continue
            
            RNA_abundance_around_selected_tss[(position_shift-tss_start_site)][1].append(RNA_coverage_site)
            
        number+=1
    
    
    if tss_strand == "-":
        #continue
        abudance_around_tss=0
    
        for position_shift in range(-100, 0):
        
            if tss_start_site+position_shift <=0 or tss_start_site+position_shift > 4471710: continue
            abudance_around_tss+=int(genome[(tss_start_site+position_shift)][3][6])
        
        
        if  abudance_around_tss <3000: continue
        if  abudance_around_tss >5000000: continue
        
        abudance_around_tss=abudance_around_tss/100
        
        check_end = tss_start_site-1000
        
        
        for position_shift in range((tss_start_site-100), tss_start_site-1000, -1):
    
            if genome[position_shift][2]>50:
        
                upstream_tss_coverage, downstream_tss_coverage=(0,1)
        
                for position_shift_tss_check in range (0, 20):
            
                    upstream_tss_coverage +=int(genome[position_shift+position_shift_tss_check][3][6])
                    downstream_tss_coverage +=int(genome[position_shift-position_shift_tss_check-1][3][6])
                if downstream_tss_coverage> upstream_tss_coverage*2 :
                    if position_shift > check_end:
                
                        check_end=position_shift+10
                        break
                        
        for position_shift in range((tss_start_site), (check_end+1), -1):
            RNA_coverage_site= int(genome[position_shift][3][6])
            #if RNA_coverage_site >5000: continue
            RNA_abundance_around_selected_tss[(tss_start_site-position_shift)][1].append(RNA_coverage_site)
            
        number+=1 
        
        
        
        
            
all_RNA_abundance_around_selected_tss=[]

for site in range(0, 1000):
    total_RNA_abundance=0
    
    if len(RNA_abundance_around_selected_tss[site][1])>0: total_RNA_abundance=sum(RNA_abundance_around_selected_tss[site][1])
    all_RNA_abundance_around_selected_tss.append([site,total_RNA_abundance ])
    
all_RNA_abundance_array=numpy.array(all_RNA_abundance_around_selected_tss)


x_avg, y_avg= numpy.array(all_RNA_abundance_array).T

#print(x_avg, y_avg)
plt.plot(x_avg, y_avg,  color="red", label="All RNA abudance", linewidth=1.5)
#plt.fill_between(x_3000, y_3000-error_3000,y_3000+error_3000, color="coral", alpha=0.25 )
#plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)    

        
        
plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)  





    
#plt.ylim((0,100))   
plt.xlim((0,1000))
plt.ylabel('Total RNA abundance')
plt.xlabel( "Distance to the operon start sites (bp)")
    #plt.xlabel( 'genome position of ' + operon_name +"(" + operon_direction + ")")
plt.savefig(input_file_name_s+"_RNA_abundance_around_selected_TSS"+"_N_"+str(number)+"_plot04122022"+".pdf", dpi=600)
    
        
print("total counted transcritps is ", number )        
                
            
            
        
    
    
    



















exit()







all_tss_3_end={}
for site in range(0, 500):
    all_tss_3_end[site]=[[],[]]
    
number=1
for genome_position in range (0,4471711):
    
    if genome[genome_position][1]>20:
        
        adjacent_tss_p=0
        adjacent_tss_n=0
        
        for position_shift in range(genome_position+1,genome_position+300):
            
            if genome[position_shift][1]>20:adjacent_tss_p+=1
            if genome[position_shift][2]>20: adjacent_tss_n+=1
            
        if adjacent_tss_p>0: continue
        
        
        
        tss_initial_expression=0
        
        
        for position_shift in range(genome_position-10,genome_position+10):
            
            tss_initial_expression+=int(int(genome[int(position_shift)][3][1]))
            
            
        if tss_initial_expression<50: continue
        
        
        
        
        
        
        Y_3end_intensity=[]
        
        for position_shift in range(genome_position,genome_position+500):
            end_intensity_p= int(int(genome[int(position_shift)][3][2]))
            if end_intensity_p >300: continue
            Y_3end_intensity.append(([(int(position_shift-genome_position)), end_intensity_p]))
            
            all_tss_3_end[(position_shift-genome_position)][1].append(end_intensity_p)
            
            
        Y_3end_intensity_array=numpy.array(Y_3end_intensity)    
        
        
        
        x_G1_p, y_G1_p = numpy.array(Y_3end_intensity_array).T
        
        number+=1
        
        plt.close()
        
        fig, axe = plt.subplots( dpi=300)
        
        #plt.plot(x_G1_p, y_G1_p,  color="red",  )
        #if operon_direction=="-": plt.plot(x_G1_n, y_G1_n,  color="blue",  )
        
            
        plt.ylabel('Intensity')
        plt.xlabel( 'Genome position of ' + str(genome_position) +"(" + "+" + ")")
        #plt.savefig('./end_output_S01/'+str(genome_position)+"_RNA_end"+".pdf")
    if genome[genome_position][2]>20:
        
        adjacent_tss_p=0
        adjacent_tss_n=0
        
        for position_shift in range(genome_position-1,genome_position-300, -1):
            if position_shift< 1: continue
            if genome[position_shift][1]>20:adjacent_tss_p+=1
            if genome[position_shift][2]>20: adjacent_tss_n+=1
            
        if adjacent_tss_n>0: continue
        
        
        
        tss_initial_expression=0
        
        
        for position_shift in range(genome_position-10,genome_position+10):
            
            tss_initial_expression+=int(int(genome[int(position_shift)][3][5]))
            
            
        if tss_initial_expression<20: continue
        
        
        
        
        
        
        Y_3end_intensity_n=[]
        
        for position_shift in range(genome_position,genome_position-500, -1):
            if position_shift< 1: continue
            end_intensity_n= int(int(genome[int(position_shift)][3][4]))
            
            Y_3end_intensity_n.append(([(int(genome_position-position_shift)), end_intensity_n]))
            if end_intensity_n >300: continue
            all_tss_3_end[(genome_position-position_shift)][1].append(end_intensity_n)
            
            
        Y_3end_intensity_n_array=numpy.array(Y_3end_intensity_n)    
        
        
        
        x_G1_n, y_G1_n = numpy.array(Y_3end_intensity_n_array).T
        
        number+=1
        
        plt.close()
        
        fig, axe = plt.subplots( dpi=300)
        
        #plt.plot(x_G1_n, y_G1_n,  color="red",  )
        #if operon_direction=="-": plt.plot(x_G1_n, y_G1_n,  color="blue",  )
        
            
        plt.ylabel('Intensity')
        plt.xlabel( 'Genome position of ' + str(genome_position) +"(" + "+" + ")")
        #plt.savefig('./end_output_S01/'+str(genome_position)+"_RNA_end"+".pdf")
        
        
        
all_end_around_TSS=[]

for site in range(0, 500):
    total_end_intensity=0
    
    if len(all_tss_3_end[site][1])>0: total_end_intensity=sum(all_tss_3_end[site][1])
    all_end_around_TSS.append([site,total_end_intensity ])
    
all_end_around_TSS_array=numpy.array(all_end_around_TSS)


x_avg, y_avg= numpy.array(all_end_around_TSS).T

#print(x_avg, y_avg)
plt.plot(x_avg, y_avg,  color="red", label="Total", linewidth=1.5)
#plt.fill_between(x_3000, y_3000-error_3000,y_3000+error_3000, color="coral", alpha=0.25 )
#plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)    

        
        
plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)  





    
#plt.ylim((0,100))   
plt.xlim((0,500))
plt.ylabel('Total intensity')
plt.xlabel( "Distance to the TSS (nt)")
    #plt.xlabel( 'genome position of ' + operon_name +"(" + operon_direction + ")")
plt.savefig(input_file_name_s+"_3ends_around_operon_TSS"+"_N_"+str(number)+"_all_TSS_plot08022021"+".pdf", dpi=600)
    
        
print("total counted transcritps is ", number )        
    








