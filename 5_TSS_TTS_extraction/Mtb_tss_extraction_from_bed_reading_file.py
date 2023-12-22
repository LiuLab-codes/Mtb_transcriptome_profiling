#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import re
import glob, os

# The inputfile should be summed *_bed_reading.txt of primary RNA SEnd-seq samples. 




def get_DNA_sequence (start_site, end_site, direction):
    #start_site -= 1
    
    sequence_output = str(sequence[start_site:end_site]) 
    if direction =="-":
        nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        sequence_output ="".join(nn[n] for n in reversed(sequence_output))
    return  sequence_output
    
genome_input_file = open ('NC_018143_TB_o.fna', 'r')
for line in  genome_input_file:
    line = line.strip() #去除前后空格
    if re.search(r"^>", line):
        continue
    sequence = line
genome_input_file.close()     



cwd = os.getcwd()
home_directory = os.getenv("HOME")
genome = {}


for num in range (0,4471711):

    genome[num] = ["","","",[],[]]  
    
print ('finish the first step of data initialization\n')





input_file_name = str(sys.argv[1])
pattern = input_file_name.split("_bed_reading")
input_file_name_s=pattern[0]
print(input_file_name_s)

bed_read_file_1 = open(input_file_name, 'r')
bed_file_coverage=0
for line in  bed_read_file_1:
    #line = line.strip() #去除前后空格
    #  1	10	0	10	0	0	9
    position_informaton = line.strip().split("\t")
    if position_informaton[0].isdigit():
        #if int(position_informaton[0]) > 210000: break
        genome[int(position_informaton[0])][2] = position_informaton
        if   1471900 < int(position_informaton[0]) < 1477050 : continue    #skip the highly expressed rRNA region.
        if int(position_informaton[0])  > 4411306 : continue
        
        
        
bed_read_file_1.close()  

print ("finish the step of reading expression data\n")

tss_extraction={}
previous_TSS=0
for genome_position in range (0, 4411306):
       # print (genome_position, vars()[sub_genome][int(genome_position)][1])
        if int(genome[genome_position][2][1]) >=10:
            if genome_position-previous_TSS <5: continue
            peak=genome_position
            for number_shift in range (0,11):
                if int(genome[genome_position+number_shift][2][1]) > int(genome[peak][2][1]):
                    peak= genome_position+number_shift
            
            total_tss_intensity=0
            for number_shift in range (-5,5):
                total_tss_intensity += int(genome[peak+number_shift][2][1])
            if total_tss_intensity <10 : continue                               #subject to change based on the overall sequencing depth
            upstream_tss_coverage, downstream_tss_coverage = (1,1)
            for number_shift in range (0,10):
                upstream_tss_coverage += int(genome[peak-number_shift-1][2][3])
                downstream_tss_coverage+= int(genome[peak+number_shift][2][3])
            
            
            if upstream_tss_coverage*1.5 > downstream_tss_coverage  : continue
            fold_change= round(downstream_tss_coverage/upstream_tss_coverage, 2)
            distance_to_previous_positive_tss= peak-previous_TSS
            previous_TSS=peak+5
            
            
            
            print(genome_position, peak,"+",total_tss_intensity, upstream_tss_coverage, downstream_tss_coverage, fold_change)
            
            peak_nucletide_1=get_DNA_sequence(peak-1,peak, "+")
            peak_nucletide_2=get_DNA_sequence(peak,peak+1, "+")
            tss_extraction[str(genome_position)+"+"]=[ peak,"+",total_tss_intensity, distance_to_previous_positive_tss, fold_change,  upstream_tss_coverage, downstream_tss_coverage, peak_nucletide_1, peak_nucletide_2, input_file_name_s]
            

previous_TSS_n= 4411306        
            
for genome_position in range (4411306, 0, -1):
       # print (genome_position, vars()[sub_genome][int(genome_position)][1])
        if int(genome[genome_position][2][5]) >=10:
            if previous_TSS_n-genome_position<5: continue
            peak=genome_position
            for number_shift in range (0,11):
                if int(genome[genome_position-number_shift][2][5]) > int(genome[peak][2][5]):
                    peak= genome_position-number_shift
            
            total_tss_intensity=0
            for number_shift in range (-5,5):
                total_tss_intensity += int(genome[peak+number_shift][2][5])
            if total_tss_intensity <10 : continue    #subject to change based on the overall sequencing depth
            upstream_tss_coverage, downstream_tss_coverage = (1,1)
            for number_shift in range (0,10):
                upstream_tss_coverage += int(genome[peak+number_shift+1][2][6])
                downstream_tss_coverage+= int(genome[peak-number_shift][2][6])
            
            
            if upstream_tss_coverage*1.5 > downstream_tss_coverage  : continue
            fold_change= round(downstream_tss_coverage/upstream_tss_coverage, 2)
            distance_to_previous_negative_tss= previous_TSS_n-peak
            previous_TSS_n=peak-5
            peak=peak
            
            
            print(genome_position, peak,"-",total_tss_intensity, upstream_tss_coverage, downstream_tss_coverage, fold_change)
            peak_nucletide_1=get_DNA_sequence(peak-1,peak, "-")
            peak_nucletide_2=get_DNA_sequence(peak-2,peak-1, "-")
            tss_extraction[str(genome_position)+"-"]=[ peak,"-",total_tss_intensity,distance_to_previous_negative_tss,fold_change,  upstream_tss_coverage, downstream_tss_coverage,peak_nucletide_1, peak_nucletide_2, input_file_name_s ]
            
            
            
print(len(tss_extraction))
            
tss_sort_position = sorted(tss_extraction.items(), key=lambda x: x[1])
output_file = open( input_file_name_s + '_tss_extraction_single_file_result.txt', 'w')
newline= "\t".join(str(i) for i in ( "genome_position", "direction", "TSS_intensity", "distance_to_previous_same_direction_tss","fold_change", "upstream_tss_coverage", "downstream_tss_coverage", "source"))
output_file.write(newline + "\n")
for tss in tss_sort_position:
    newline= "\t".join(str(i) for i in (tss[1]))
    output_file.write(newline + "\n")
    
output_file.close()
    
    
    



            
            
            
        



