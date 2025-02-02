set -ue

# To convert the merged fasta sequence data to paired end data, this program will remove the 5'end adator at last
# You need to install the cutadapt, FASTX-Toolkit and bowtie2 software package to run this script.

# get name inforamtion

fa_file=$1
name=$2

#remove the adaptor duplicated reads.
# remove the reads containing the duplicated adaptors.
# This script is designed to split the fasta into paired-end SEnd-seq reads. The 3' adaptor required is "/5Phos/rNrNrNrNrArArCrCrUrGrCrUrArUrCrArArCrUrG/3ddC/". And the sequence information here is subjected to change based on different 3' RNA adaptors.


grep -vE "(ATAGCAGGTT).*(ATAGCAGGTT)|(AACCTGCTAT).*(AACCTGCTAT)" $fa_file | grep -vE "(ACGTGGAGAG).*(ACGTGGAGAG)|(CTCTCCACGT).*(CTCTCCACGT)" | grep -vE "GGGCAGTTGATAGCAGGTTGCTAT|ATAGCAACCTGCTATCAACTGCCC" | grep -vE "GGGCAGTTGATAGCAGGAAGC|GCTTCCTGCTATCAACTGCCC" |grep -E "GGGCAGTTGATAGCAGGTT|AACCTGCTATCAACTGCCC" -B 1 | grep -v "^--$" > $name.3L.fa

#seperate the forward adaptor and reverse adaptor, and then reverse complementary change the reverse adaptor congtaing reads. and then combine the two files togehter.


grep -E "GGGCAGTTGATAGCAGGTT" -B 1 $name.3L.fa | grep -v "^--$" > $name.f.fa
grep -vE "GGGCAGTTGATAGCAGGTT"  $name.3L.fa | grep -E "AACCTGCTATCAACTGCCC" -B 1 | grep -v "^--$" > $name.r.fa
fastx_reverse_complement -i $name.r.fa -o $name.r_RC.fa

cat $name.f.fa $name.r_RC.fa > $name.c2.fa 
rm $name.f.fa $name.r_RC.fa $name.r.fa

# seperate the adaptor containing reads into two parts, one is for the 5' end and one is for the 3'end

sed "s/GGGCAGTTGATAGCAGGTT/XXXXX\\`echo -e '\n\r'`YYYYY/g"  $name.c2.fa  > $name.sep.c2.fa

# seperate the 5' end into another file, and remove unnecessary parts

grep -E "XXXXX" -B 1  $name.sep.c2.fa  |sed "s/XXXXX//g"  | grep -v '^--$' > $name.sep.c2.P1.fa

# seperate the 3' end parts into another file.

grep -vE "XXXXX"  $name.sep.c2.fa | grep -E "YYYYY" -B1 | sed "s/YYYYY//g" |grep -v '^--$' > $name.sep.c2.P2.fa

dos2unix $name.sep.c2.P2.fa

# use the cutadapt software to remove possible 5' adaptor and remove too short reads

cutadapt  -e 0.1    -U 4 -m 15  --trim-n -o $name.sep.c2.P1_cut.fa -p $name.sep.c2.P2_cut.fa $name.sep.c2.P1.fa  $name.sep.c2.P2.fa

# cut the final reads to less than 75bp and then make the sequence of reads to starndard data format.


cutadapt   -o $name.sep.c2.P2_cut_75bp.fa $name.sep.c2.P2_cut.fa

fastx_reverse_complement -i $name.sep.c2.P2_cut_75bp.fa -o $name.sep.c2.P2_cut_75bp_RC.fa
fastx_reverse_complement -i $name.sep.c2.P1_cut.fa -o $name.sep.c2.P1_cut_RC.fa
cutadapt   -o $name.sep.c2.P1_cut_RC_75bp.fa $name.sep.c2.P1_cut_RC.fa


# save the final files to the fold of final.

file="final"
if [ -d "$file" ]
then
	echo "$file exists and the new files will be added into final."
else
	mkdir final
fi

 cp *RC_75bp.fa  *75bp_RC.fa final/
 rm $name.c2.fa  $name.sep.c2.fa $name.sep.c2.P1.fa $name.sep.c2.P2.fa $name.sep.c2.P1_cut.fa $name.sep.c2.P2_cut.fa  $name.sep.c2.P2_cut_75bp.fa $name.sep.c2.P1_cut_RC.fa $name.sep.c2.P1_cut_RC_75bp.fa   $name.sep.c2.P2_cut_75bp_RC.fa

 bowtie2 -p8 --no-mixed  -f --very-sensitive-local   -X 10000 --ff -x ~/ref_E.coli/NC_000913_spikeRNA -1 ./final/$name.sep.c2.P1_cut_RC_75bp.fa  -2 ./final/$name.sep.c2.P2_cut_75bp_RC.fa  -S ./final/$name\_3L.bowtie2mapping.sam 
