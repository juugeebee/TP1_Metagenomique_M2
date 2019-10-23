#!/bin/bash

#./16S_workflow.sh fastq/ fastqc_out/

dirfq=$1
outdir=$2

mkdir -p $2/fastqc_avanttrim
mkdir -p $2/fastqc_aprestrim
mkdir -p fasta

####deziper les fastq
gunzip $1/*.gz

####fastqc
fastqc $1/*_R1.fastq $1/*_R2.fastq -o $2/fastqc_avanttrim
for i in $list_fastq; do
    fastqc -o $2/fastqc_avanttrim --noextract fastq/$i
done

for file in $(ls $1/*_R1.fastq)
do
	nameR1="$file"
	nameR2=$(echo $file | sed "s/R1/R2/g")
	samplename=$(echo $file | cut -d/ -f 2 | cut -d_ -f 1)
	
	#trimm
	java -jar soft/AlienTrimmer.jar -if $nameR1 -ir $nameR2 -c databases/contaminants.fasta -q 20
	
	#fastqc
	fastqc $1/*_R1.fastq.at.fq $1/*_R2.fastq.at.fq -o $2/fastqc_aprestrim
	
	#Pairer les séquence F et R > Merger
	vsearch --fastq_mergepairs $1/*_R1.fastq.at.fq --reverse $1/*_R2.fastq.at.fq --fastaout fasta/$samplename.fasta --label_suffix "$s;sample=$samplename;" --fastq_minovlen 40 --fastq_maxdiffs 15   
	
done





