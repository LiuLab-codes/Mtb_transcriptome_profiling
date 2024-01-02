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



gene_cover_start_left = 3568560     #specify the start site of the targeting region.
gene_cover_start_right = 3569050    #specify the end site of the targeting region.
target_direction = "-"              #specify the direction of the targeting region.
length= gene_cover_start_right -  gene_cover_start_left
TSS_transcript_infor = {}
genome = {}
genome_seq_signal={}
cwd = os.getcwd()

home_directory = os.getenv("HOME")
for num in range (0,4471711):

    genome[num] = [0,0,0,"","",0,0,0]
    
print ('finish the first step of data initialization\n')


gff_input_file = open('NC_018143.gff', 'r')

gff_infor = {}
for line_gff in  gff_input_file:
    #line_bed = line_bed.strip() #去除前后空格
    position_informaton_gff = line_gff.strip().split("\t")
    pattern = "name=(.*?);"
    gene_name=re.search(pattern, position_informaton_gff[8]).group(1)
    
    gene_start_site, gene_end_site, direction = (int(position_informaton_gff[3]), int(position_informaton_gff[4]), position_informaton_gff[6])
    
    
    for  position in range (gene_start_site, gene_end_site+1):
        genome[position][0]=gene_name
        
    #print(gene_name, gene_start_site, gene_end_site, direction)
    if gene_name in gff_infor : print ("wrong gene name as a key")
    gff_infor[gene_name] = [gene_name, gene_start_site, gene_end_site, direction]
    
    
gff_input_file.close()














input_file_name_1 = str(sys.argv[1])
pattern = input_file_name_1.split("_bed_reading")
input_file_name_s_1=pattern[0]

bed_read_file_1 = open(input_file_name_1, 'r')


for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        genome[int(position_informaton[0])][3] = position_informaton

        
bed_read_file_1.close()   






X_nucletide_posiiton = []

Y_1_coverage_p = []
Y_1_coverage_n = []

G1_coverage_origin_p=[]
G1_coverage_origin_n=[]

y_max_p, y_max_n=(0,0)
print(gene_cover_start_left,gene_cover_start_right)
for screen_position in range (gene_cover_start_left, (gene_cover_start_right+1)):
        G1_coverage_origin_p.append([(int(screen_position)), int(genome[int(screen_position)][3][3])])
        G1_coverage_origin_n.append([(int(screen_position)), int(genome[int(screen_position)][3][6])])
        if int(genome[int(screen_position)][3][3])>y_max_p: y_max_p=int(genome[int(screen_position)][3][3])
        if int(genome[int(screen_position)][3][6]) > y_max_n: y_max_n=int(genome[int(screen_position)][3][6])



G1_coverage_p_array = numpy.array(G1_coverage_origin_p)
G1_coverage_n_array = numpy.array(G1_coverage_origin_n)

x_G1_p, y_G1_p = numpy.array(G1_coverage_p_array).T
x_G1_n, y_G1_n = numpy.array(G1_coverage_n_array).T


plt.close()
fig, axe = plt.subplots( dpi=300)

coding_gene=[]

for position in range(gene_cover_start_left,gene_cover_start_right+1) :
     
    gene_name =genome[position][0]
    
    if gene_name == 0 or gene_name in coding_gene : continue
    coding_gene.append(gene_name)
    
    
    

gene_serial=1
if len(coding_gene)>=1: 
    for gene in coding_gene:
    
        gene_start_site, gene_end_site, gene_direction = (int(gff_infor[gene][1]), int(gff_infor[gene][2]),gff_infor[gene][3])
        
        if target_direction == "+": 
            x_gene, y_gene = numpy.array([[(int(gene_start_site)), -4*gene_serial], [int(gene_end_site)-length/30, -4*gene_serial]]).T
            color_gene = "coral"
            gene_line=plt.plot(x_gene, y_gene,  color=color_gene,  linewidth=2)

            axe.arrow(int(gene_end_site)-length/30, -4*gene_serial, length/30*0.2, 0,linewidth=2, head_width=y_max_p/90, head_length=length/30*0.6, color=color_gene)
            axe.text(gene_start_site+int((gene_end_site-gene_start_site)/3), -4*gene_serial+y_max_p/80, str(gene), fontsize=7)
        elif target_direction == "-": 
            
            
            x_gene, y_gene = numpy.array([[(int(gene_start_site))+length/30, -4*gene_serial], [int(gene_end_site), -4*gene_serial]]).T
            color_gene = "royalblue"
            gene_line=plt.plot(x_gene, y_gene,  color=color_gene,  linewidth=2)

            axe.arrow(int(gene_start_site)+length/30, -4*gene_serial, -length/30*0.2, 0,linewidth=2, head_width=y_max_n/90, head_length=length/30*0.6,  color=color_gene)
            axe.text(gene_start_site+int((gene_end_site-gene_start_site)/3), -4*gene_serial+y_max_n/80, str(gene), fontsize=7)
            
            
        #x_gene, y_gene = numpy.array([[(int(gene_start_site)), -4*gene_serial], [int(gene_end_site), -4*gene_serial]]).T
        #color_gene = "orange"
        #if gene_direction == "-": color_gene = "cyan"
        #plt.plot(x_gene, y_gene,  color=color_gene,  linewidth=2, alpha=0.5)
    
        #axe.text(int((gene_start_site+gene_end_site)/2), -4*gene_serial, str(gene), fontsize=8)
        gene_serial+=1
    







if target_direction=="+":   plt.plot(x_G1_p, y_G1_p,  color="red",  )
if target_direction=="-": plt.plot(x_G1_n, y_G1_n,  color="blue",  )






plt.ylabel('Coverage')
plt.xlabel( 'Genome position')
plt.savefig(str(gene_cover_start_left)+"_RNA_coverage"+".pdf", dpi=300)


exit()
    
    
    








