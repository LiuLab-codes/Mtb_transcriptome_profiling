#!/usr/bin/env perl -w 
use strict; 

use File::HomeDir;
my $fileSpec = File::HomeDir->my_home;


#Requirements: Before running this script, ensure you have installed Trim Galore (https://github.com/FelixKrueger/TrimGalore), FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), Cutadapt (https://cutadapt.readthedocs.io/en/stable/index.html), and FASTX-Toolkit (http://hannonlab.cshl.edu/fastx_toolkit/). Additionally, specify the directory path of Cutadapt in line 26 of the script. Note: This script is not compatible with Apple silicon chips.

#The phrase "R1_001" in lines 20 and 21 should be modified according to the naming convention of your raw data files. In this example, we are using files ending with "_R1_001.fastq" or "_R2_001.fastq".
#Place this script in the same directory as your raw data. Execute it by running the command 'perl pair_end_fastq_merge.pl' in your command line terminal." The file with the extension '*.n.fa' is anticipated to be the output file.
# ß
my $directory = `pwd`;
chomp $directory;
#print "$directory";
opendir(DIR, $directory) or die $!;

my $file_name_1;
while (my $file = readdir(DIR)) {
  next if  $file =~ /\.$/;
 if ($file =~ /R1_001/){
 	$file=~ /_R1_001/;
	 $file_name_1 = $`;
	print "$file_name_1\n";
	#&pair_end_fq_check("$file_name_1\_R1_001.fastq","$file_name_1\_R2_001.fastq"); # Here, should you experience any formatting errors with your paired-end (PE) data, you can activate this command to perform a check on your input data files
	system("trim_galore --path_to_cutadapt  /Users/Xju/opt/miniconda3/bin/cutadapt -q 20  --length 30 --paired --retain_unpaired $file_name_1\_R1_001.fastq $file_name_1\_R2_001.fastq ");
	system("flash -r 150 -f 250 -s 150 $file_name_1\_R1_001_val_1.fq  $file_name_1\_R2_001_val_2.fq  -o ./$file_name_1 2>&1 | tee $file_name_1\_flash.log");
	system(" cat $file_name_1.extendedFrags.fastq $file_name_1.notCombined_1.fastq $file_name_1.notCombined_2.fastq $file_name_1\_R1_001_unpaired_1.fq $file_name_1\_R2_001_unpaired_2.fq > $file_name_1.fq ");
	system (" fastq_to_fasta -Q33 -i  $file_name_1.fq -o $file_name_1.fa ");
	
	system (" fastx_renamer -n COUNT -i $file_name_1.fa -o $file_name_1.n.fa");
	system (" rm $file_name_1\_R1_001.fastq $file_name_1\_R2_001.fastq $file_name_1\_R1_001_val_1.fq  $file_name_1\_R2_001_val_2.fq  $file_name_1.extendedFrags.fastq $file_name_1.notCombined_1.fastq $file_name_1.notCombined_2.fastq $file_name_1\_R1_001_unpaired_1.fq $file_name_1\_R2_001_unpaired_2.fq  $file_name_1.fa  $file_name_1.fq ");
 };
  
 

}





#&pair_end_fq_check("total-RD-TB-WT-log_S15_L001_R1_001.fastq","total-RD-TB-WT-log_S15_L003_R2_001.fastq");








sub pair_end_fq_check {
my ($paired_fq_file_1, $paired_fq_file_2) = @_; #Open bed file and sam file
print "$paired_fq_file_1\t$paired_fq_file_2\n";

my %fasta_name;





open (fastq_file_1, "./$paired_fq_file_1") or die "can't open the input fastq_file\n";
my $previous_line_number =-3;
my $line_number =0;
my @lines_file_1;
my $target_line =4;
while (<fastq_file_1>){
	$line_number +=1;
	$target_line +=1;
chomp;

my $line_1 = "$_";
if (($line_number == ($previous_line_number+4)) && ($line_1 =~ /^\Q@\E/)){
	my @line_name= split(/ /);
	#print "$line_name[0]\n";
	$fasta_name{$line_number}=$line_name[0];
		$target_line =1;
		
	$previous_line_number +=4;
}




}

close fastq_file_1;


open (fastq_file_2, "./$paired_fq_file_2") or die "can't open the input fastq_file\n";
 $previous_line_number =-3;
 $line_number =0;
my @lines_file_2;
$target_line =4;
while (<fastq_file_2>){
	$line_number +=1;
	$target_line +=1;
chomp;

my $line_2 = "$_";
if (($line_number == ($previous_line_number+4)) && ($line_2 =~ /^\Q@\E/)){
	my @line_name= split(/ /);
	#print "$line_name[0]\n";
	#print "$fasta_name{$line_number}\t$line_name[0]\n";
	print "$paired_fq_file_2 is wrong\n" unless ($fasta_name{$line_number} eq $line_name[0]);
	last unless $fasta_name{$line_number} eq $line_name[0];	#print "$line_name[0]\n";
	}
	$previous_line_number +=4;
}






close fastq_file_2;
}
