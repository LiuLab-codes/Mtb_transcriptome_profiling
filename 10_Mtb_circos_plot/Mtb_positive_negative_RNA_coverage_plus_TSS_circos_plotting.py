
#!/use/bin/env python
# -*- coding: UTF-8 -*-
#%matplotlib inline
import sys 
import os.path
import glob, os
import numpy as np
sys.path.append("./")
from pycircos import *
from Bio import SeqIO
import time
import math

genome = {}
if __name__ == "__main__":
    
    

    #input_file_name = os.path.splitext(str(sys.argv[1]))[0]
    home_directory = os.getenv("HOME")
    for num in range (0,4411708):
        genome[num] = ["",0,0,"","",0,0,0]
        #genome_seq_signal[num] = [num,0,0,0,0,0,0]
        #if num > 220000 : break

    print ('finish the first step of data initialization\n')

    
            
    
    bed_read_file_2 =  open(str(sys.argv[1]), 'r')

    for line in  bed_read_file_2:  #total RNA data
        
    
        position_informaton = line.strip().split("\t")
        if position_informaton[0].isdigit():
            if int(position_informaton[0]) <4411708:
                genome[int(position_informaton[0])][4] = position_informaton

        
    bed_read_file_2.close()
    
    TSS_input_file= open ('Mtb_TSS_final.txt', 'r')
    for line_tss in  TSS_input_file:
        position_informaton_tss = line_tss.strip().split("\t")
        #print(position_informaton_tss)
        if len(position_informaton_tss)> 2 and position_informaton_tss[0].isdigit():
            position, direction, expression_level = (int(position_informaton_tss[0]),position_informaton_tss[1], int(position_informaton_tss[2]))
            sum_around_positive_tss, sum_around_negative_tss = (0,0)
            
            if direction =="+":
                for num_shift in range(-5,5):
                    if position+num_shift >= 4411708: continue
                    sum_around_positive_tss+=int(genome[position+num_shift][4][1])
                genome[position][5]=sum_around_positive_tss
                print(genome[position][5])
                
            if direction =="-":
                for num_shift in range(-5,5):
                    if position+num_shift >= 4411708: continue
                    sum_around_negative_tss+=int(genome[position+num_shift][4][5])
                genome[position][6]=sum_around_negative_tss
                print(genome[position][6])
    TSS_input_file.close()
    
    

    print(genome[4411254])
    #print ("finish the step of loading the bed_reading information\n")


    positive_coverages = []
    negative_coverages = []
    TSSs_positive= []
    TSSs_negative = []

    window_size=1
    slide_size=None
    if slide_size is None:
        slide_size = window_size
    for i in range(0, 4411708, slide_size):
    
        positive_coverage=0
        negative_coverage=0
        TSS_positive =0
        TSS_negative =0
    
        if (i+slide_size) > 4411700 :
            for position in range (i, 4411708):
                positive_coverage += int(genome[position][4][3])
                negative_coverage += int(genome[position][4][6])
                TSS_positive += int(genome[position][5])
                TSS_negative += int(genome[position][6])
            positive_coverage = round((positive_coverage/(4411709-i)),2)
            negative_coverage = round((negative_coverage/(4411709-i)), 2)
            
            
            
            
            
            
        else:
        
            for position in range (i, (i+slide_size)):
                positive_coverage += int(genome[position][4][3])
                negative_coverage += int(genome[position][4][6])
                TSS_positive += int(genome[position][5])
                TSS_negative += int(genome[position][6])
            #positive_coverage = round((positive_coverage/slide_size),2)
            #negative_coverage = round((negative_coverage/slide_size), 2)
        
    
        if positive_coverage <=1 : 
            positive_coverage=0
        else:
            positive_coverage = round(math.log10(positive_coverage),3)
        #print(i, positive_coverage)

        positive_coverages.append(positive_coverage)
    
    
        if negative_coverage <=1 : 
            negative_coverage=0
        else:
            negative_coverage = round(math.log10(negative_coverage),3)
        #print(i, positive_coverage)

        negative_coverages.append(negative_coverage)
        
        
        if TSS_positive <=0 :
            TSS_positive =0 
        else:
            TSS_positive = round(math.log10(TSS_positive),3)
            
            
        if TSS_negative <=0 :
            TSS_negative =0 
        else:
            TSS_negative = round(math.log10(TSS_negative),3)
        
            
            
        
        
        
        TSSs_positive.append(TSS_positive)
        TSSs_negative.append(TSS_negative)
            
            
            
            
            
            
            
            
    positive_coverages = np.array(positive_coverages)
    negative_coverages = np.array(negative_coverages)
    TSSs_positive = np.array(TSSs_positive)
    TSSs_negative = np.array(TSSs_negative)
    
    print ("coverage and TSS max value ")
    print(np.max(positive_coverages))
    print(np.max(negative_coverages))
    print(np.max(TSSs_positive))
    print(np.max(TSSs_negative))
    print ("coverage and TSS mean value ")
    print(np.mean(positive_coverages))
    print(np.mean(negative_coverages))
    print(np.mean(TSSs_positive))
    print(np.mean(TSSs_negative))
    print ("coverage and TSS min value ")
    print(np.min(positive_coverages))
    print(np.min(negative_coverages))
    print(np.min(TSSs_positive))
    print(np.min(TSSs_negative)) 
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    gbk = SeqIO.parse("GCF_000277735.2_ASM27773v2_genomic.gbff","genbank")
    #print (gbk)
    
    
    
    
    
    
    
    
    #exit()
    gcircle = Gcircle()
    gcircle.interspace = 0.0
    gcircle.read_locus(gbk, bottom=900, height=0, linewidth=0) 
    gcircle.set_locus()
    locus_names = list(gcircle.locus_dict.keys())
    #print(locus_names)
    
    
    
    
    #positive coverage of SEnd-seq
    
    gcircle.line_plot('NC_018143.2', positive_coverages, bottom=700, height=150, facecolor="r", linewidth=0, fill=True) 
    gcircle.line_plot('NC_018143.2', negative_coverages, bottom=700, height=-150, facecolor="b", linewidth=0, fill=True)
    #gcircle.bar_plot('NC_018143.2', positive_coverages, bottom=600, height=100, facecolor="#888888")
    
    
    
    #TSS of SEND-seq
    gcircle.line_plot('NC_018143.2', TSSs_positive, bottom=500, height=100, facecolor="r", linewidth=0.1, edgecolor="r",fill=True) 
    gcircle.line_plot('NC_018143.2', TSSs_negative, bottom=500, height=-100, facecolor="b", linewidth=0.1, edgecolor="b",fill=True)
    
    
    #GC skews
    
    #Visualization of gc skew by fill plot
    #gc_skew = gcircle.calc_gcskew('NC_018143.2', window_size=1000)
    #gc_skew_positive = np.copy(gc_skew) 
    #gc_skew_positive[gc_skew_positive<0] = 0 #Get positive value of gc skew
    #gc_skew_negative = np.copy(gc_skew)
    #gc_skew_negative[gc_skew_positive>0] = 0 #Get negative value of gc skew
    #gcircle.line_plot('NC_018143.2', gc_skew_positive, bottom=500, height=100, facecolor="r", linewidth=0, fill=True) #plot positive gc skew
    #gcircle.line_plot('NC_018143.2', -1.0 * gc_skew_negative, bottom=500, height=-100, facecolor="b", linewidth=0, fill=True) #plot negative gc skew
    gcircle.set_spine('NC_018143.2', 500, linewidth=0.01) #set spine
    gcircle.set_spine('NC_018143.2', 700, linewidth=0.01) #set spine
    #fig = gcircle.save()
    
    
    #Visualization of gc ratio by barplot
    
    #gc_ratio = gcircle.calc_gcratio('NC_018143.2', window_size=1000)
    #gcircle.bar_plot('NC_018143.2', gc_ratio, bottom=200, height=100, facecolor="#888888")
    
    
    #print(gc_ratio)
    #Visualization of CDS 
    
    
    
    
    
    
    
    
    
    gcircle.plot_features('NC_018143.2', bottom=840, height=20, facecolor= "#6388b4" ,requirement=lambda x:x.location.strand==-1)
    gcircle.plot_features('NC_018143.2', bottom=860, height=20, facecolor="#ef6f6a", requirement=lambda x:x.location.strand==1)
    
    
    
    
    
    ticks= [1]
    

    for n in range(1,100):
        if n*200000 < 4411501:
            ticks.append(n*200000)
    
    gcircle.tick_plot('NC_018143.2', ticks, bottom=900)
    timestr = time.strftime("%Y%m%d-%H%M%S")
    print (timestr)
    gcircle.save(file_name="m.tb_"+timestr)
