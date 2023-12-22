#!/use/bin/env python
# -*- coding: UTF-8 -*-
import sys
import os.path
import re



# The inputfile should be summed *_bed_reading.txt of total RNA SEnd-seq samples. Verification of the final results with primary RNA SEnd-seq is necessary.

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
#for num in range (0,271711):
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
        if   1471900 < int(position_informaton[0]) < 1477050 : genome[int(position_informaton[0])][2]=[int(position_informaton[0]),0,0,0,0,0,0]
        if int(position_informaton[0])  > 4411306 : continue
        #bed_file_coverage += int(genome[int(position_informaton[0])][2][3])
        #bed_file_coverage += int(genome[int(position_informaton[0])][2][6])
        
        

        
bed_read_file_1.close()  

print ("finish the step of reading expression data\n")

tts_extraction={}
previous_TTS=0
for genome_position in range (0, 4411306):
       # print (genome_position, vars()[sub_genome][int(genome_position)][1])
        if int(genome[genome_position][2][2]) >10:
            if genome_position-previous_TTS <10: continue
            peak=genome_position
            for number_shift in range (0,11):
                if int(genome[genome_position+number_shift][2][2]) > int(genome[peak][2][2]):
                    peak= genome_position+number_shift
            
            total_tts_intensity=0
            for number_shift in range (-5,5):
                total_tts_intensity += int(genome[peak+number_shift][2][2])
            if total_tts_intensity <20 : continue                    
            upstream_tts_coverage, downstream_tts_coverage = (1,1)
            for number_shift in range (0,10):
                upstream_tts_coverage += int(genome[peak-number_shift][2][3])
                downstream_tts_coverage+= int(genome[peak+number_shift+1][2][3])
            
            
            #if upstream_tts_coverage < downstream_tts_coverage*1.5  : continue
            fold_change= round(downstream_tts_coverage/upstream_tts_coverage, 3)
            distance_to_previous_positive_tts= peak-previous_TTS
            previous_TTS=peak+5
            
            termination_efficiency= (1-fold_change)*100
            
            
            if termination_efficiency< 50: continue
            print(peak,"+",total_tts_intensity, upstream_tts_coverage, downstream_tts_coverage, termination_efficiency)
            tts_extraction[str(genome_position)+"+"]=[ peak,"+",total_tts_intensity, termination_efficiency,  upstream_tts_coverage, downstream_tts_coverage, input_file_name_s]
            

previous_TTS_n= 4411306        
            
for genome_position in range (4411306, 0, -1):
       # print (genome_position, vars()[sub_genome][int(genome_position)][1])
        if int(genome[genome_position][2][4]) >10:
            if previous_TTS_n-genome_position<10: continue
            peak=genome_position
            for number_shift in range (0,11):
                if int(genome[genome_position-number_shift][2][4]) > int(genome[peak][2][4]):
                    peak= genome_position+number_shift
            
            total_tts_intensity=0
            for number_shift in range (-5,5):
                total_tts_intensity += int(genome[peak+number_shift][2][4])
            if total_tts_intensity <20 : continue         #subject to change based on the overall sequencing depth
            upstream_tts_coverage, downstream_tts_coverage = (1,1)
            for number_shift in range (0,10):
                upstream_tts_coverage += int(genome[peak+number_shift][2][6])
                downstream_tts_coverage+= int(genome[peak-number_shift-1][2][6])
            
            
            #if upstream_tts_coverage < downstream_tts_coverage*1.5   : continue
            fold_change= round(downstream_tts_coverage/upstream_tts_coverage, 2)
            distance_to_previous_negative_tts= previous_TTS_n-peak
            
            termination_efficiency= (1-fold_change)*100
            
            if termination_efficiency< 50: continue
            previous_TTS_n=peak-5
            #peak=peak
            
            
            print(peak,"-",total_tts_intensity, upstream_tts_coverage, downstream_tts_coverage, termination_efficiency)
            
            tts_extraction[str(genome_position)+"-"]=[ peak,"-",total_tts_intensity,termination_efficiency,  upstream_tts_coverage, downstream_tts_coverage,input_file_name_s ]
            
            
            
print(len(tts_extraction))
            
tts_sort_position = sorted(tts_extraction.items(), key=lambda x: x[1])
output_file = open( input_file_name_s + '_tts_extraction_single_file.txt', 'w')
newline= "\t".join(str(i) for i in ( "genome_position", "direction", "TTS_intensity","Termination_efficiency", "upstream_tts_coverage", "downstream_tts_coverage", "source"))
output_file.write(newline + "\n")
for tts in tts_sort_position:
    newline= "\t".join(str(i) for i in (tts[1]))
    output_file.write(newline + "\n")
    
output_file.close()
    
    
    



            
            
            
        



