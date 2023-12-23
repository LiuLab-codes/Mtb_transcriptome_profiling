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
import statistics 
from statistics import mode 



genome = {}

cwd = os.getcwd()


for num in range (0,4471711):

    genome[num] = ["",0,0,"","","",[],[]]
    
    
print ('finish the first step of data initialization\n')



input_file_name = str(sys.argv[1])
pattern = input_file_name.split(".sorted_bed_reading")
input_file_name_s=pattern[0]

spike_in_RNA_coverage=1
RNA_coverage=0
bed_read_file_1 = open(input_file_name, 'r')
for line in  bed_read_file_1:
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        
        genome[int(position_informaton[0])][3] = position_informaton
    
    
        
bed_read_file_1.close()   



average_1000_RNA={}

for site in range(-100, 701):
    average_1000_RNA[site]=[[],[]]
   



    
number=0

tss_select_file  = open("select_tss_without_downstream_TSSs_in_short_region.txt", 'r')

for tss_line in tss_select_file:
    tss_information=tss_line.strip().split("\t")
    tss_start_site, tss_strand, tss_intensity, tss_category, upstream_positive_tss_distance, upstream_negative_tss_distance, downstream_positive_tss_distance, downstream_negative_tss_distance=tss_information
    
    
    
    tss_start_site=int(tss_start_site)
    if tss_strand == "+":
        
        start_site_expression=1
        for  screen_position in range (tss_start_site, (tss_start_site+200)):
            
            start_site_expression+=int(genome[int(screen_position)][3][3])
        
        if start_site_expression<2000: continue  # The criteria for filtering are subject to change.
            
        max_covreage_RNA=1
        for  screen_position in range (tss_start_site-100, (tss_start_site+701)):
            
            
            if int(genome[int(screen_position)][3][3])>max_covreage_RNA: max_covreage_RNA= int(genome[int(screen_position)][3][3])
            
        
        
        
        for screen_position in range  (tss_start_site-100, (tss_start_site+701)):
            
            
            average_1000_RNA[(int(screen_position)-tss_start_site)][1].append(int(genome[int(screen_position)][3][3])/max_covreage_RNA*100)
        
        
        number+=1
    
    
    if tss_strand == "-":
        
        
        max_covreage_RNA=1
        
        
        start_site_expression=1
        for  screen_position in range (tss_start_site, (tss_start_site-200), -1):
            
            start_site_expression+=int(genome[int(screen_position)][3][6])
        
        if start_site_expression<2000: continue  # The criteria for filtering are subject to change.
            
        
        
        
        
        for  screen_position in range (tss_start_site+100, (tss_start_site-701), -1):
            
            
            if int(genome[int(screen_position)][3][6])>max_covreage_RNA: max_covreage_RNA= int(genome[int(screen_position)][3][6])
            
            
        X_nucletide_posiiton = []
        Y_1_coverage_n = []
        
        G1_coverage_origin_n=[]
        
        
        
        for screen_position in range  (tss_start_site+100, (tss_start_site-701), -1):
            
            
            average_1000_RNA[-(int(screen_position)-tss_start_site)][1].append(int(genome[int(screen_position)][3][6])/max_covreage_RNA*100)
        
        number+=1
    
    

all_average_1000_RNA=[]
for site in range(-100, 701):
    average_coverage_1000_RNA=0
    if len(average_1000_RNA[site][1])>0: average_coverage_1000_RNA =   numpy.median(average_1000_RNA[site][1])        
    
    error_1000_RNA=0
    if len(average_1000_RNA[site][1])>2:error_1000_RNA=statistics.stdev(average_1000_RNA[site][1])
    
    all_average_1000_RNA.append([site, average_coverage_1000_RNA, error_1000_RNA ])
    

all_average_1000_RNA_array=numpy.array(all_average_1000_RNA)

x_1000_RNA, y_1000_RNA, error_1000_RNA = numpy.array(all_average_1000_RNA_array).T

plt.plot(x_1000_RNA, y_1000_RNA,  color="purple", label="RNA", linewidth=1.5)
plt.fill_between(x_1000_RNA, y_1000_RNA-error_1000_RNA,y_1000_RNA+error_1000_RNA, color="pink", alpha=0.25 ) 


plt.legend(bbox_to_anchor=(0.7, 1), loc='upper left',  fontsize=8)  

plt.ylim((0,120))   
plt.xlim((-100,700))
plt.ylabel('Normalized intensity (%)')
plt.xlabel( "Distance to the TSSs (nt)")
   
    
plt.savefig(input_file_name_s+"_RNA_abundance_downstream_selected_tss"+"_N_"+str(number)+"_select_tss_plot"+".pdf", dpi=600)
    
