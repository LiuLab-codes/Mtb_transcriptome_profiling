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
gene_cover_start_left = 600398
gene_cover_start_right = 601854
target_direction = "+"
TSS_transcript_infor = {}
genome = {}
genome_seq_signal={}
cwd = os.getcwd()
#input_file_name = os.path.splitext(str(sys.argv[1]))[0]
home_directory = os.getenv("HOME")
for num in range (0,4471711):
#for num in range (0,271711):
    genome[num] = ["",0,0,"","",0,0,0]
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


TSS_input_file =  open('TSS_collection06212021.txt', 'r')

for line in  TSS_input_file:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        #if int(position_informaton[0]) >21000:
            #continue
    
        if (len(position_informaton) >3):
            if (str(position_informaton[1]) == "+"):
                genome[int(position_informaton[0])][1]=position_informaton[2]
                
                
            elif (str(position_informaton[1]) == "-"):
                genome[int(position_informaton[0])][2]=position_informaton[2]
        
    #print (position_informaton[0], position_informaton[1], position_informaton[2])
TSS_input_file.close()


print ("finish the step of loading the TSS information\n")

bed_read_file_1 = open('S05_total_RD_TB_WT_log_bed_reading.txt', 'r')

for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        genome[int(position_informaton[0])][3] = position_informaton

        
bed_read_file_1.close()   

gff_input_file = open('NC_018143.gff', 'r')

gff_infor = {}
for line_gff in  gff_input_file:
    #line_bed = line_bed.strip() #去除前后空格
    position_informaton_gff = line_gff.strip().split("\t")
    pattern = "name=(.*?);"
    gene_name=re.search(pattern, position_informaton_gff[8]).group(1)
    
    gene_start_site, gene_end_site, direction = (int(position_informaton_gff[3]), int(position_informaton_gff[4]), position_informaton_gff[6])
    #print(gene_name, gene_start_site, gene_end_site, direction)
    if gene_name in gff_infor : print ("wrong gene name as a key")
    gff_infor[gene_name] = [gene_name, gene_start_site, gene_end_site, direction]
    
    
gff_input_file.close()


operon_input_file = open('operon_annoation06252021.txt', 'r')

or line in  operon_input_file:
    position_informaton = line.strip().split("\t")
    operon_infor[position_informaton[0]]=position_informaton
    if not position_informaton[1].isdigit(): continue
    
    operon_seral_number, operon_start_site, operon_end_site,  operon_length, operon_direction, coding_start, coding_end, essential_check, included_gene_number, *coding_gene =position_informaton
    
    
    
    





    
    








