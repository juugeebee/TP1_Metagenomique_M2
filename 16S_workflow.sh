#!/bin/bash
chmod -x 16S_workflow.sh

folderin=$1
folderout=$2

mkdir $2
mkdir $2/fastqc_avanttrim
mkdir $2/fastqc_aprestrim
mkdir $2/trimmed
mkdir $2/fasta

#### Deziper les fastq
gunzip $1/*.gz

######## PARTIE 1 ########

#### Reportez les distributions de qualités avant filtrage qualité à l’aide de fastqc
fastqc $1/*_R1.fastq $1/*_R2.fastq -o $2/fastqc_avanttrim

for i in $list_fastq; do
    fastqc -o $2/fastqc_avanttrim --noextract fastq/$i
done


#### Executer le script d'Alientrimmer
cd soft
./JarMaker.sh AlienTrimmer.java
cd ..


for file in $(ls $1/*_R1.fastq)
do
	nameR1="$file"
	out_nameR1=$(echo $nameR1 | cut -d/ -f 2 | cut -d. -f 1)
	
	nameR2=$(echo $file | sed "s/R1/R2/g")
	out_nameR2=$(echo $nameR2 | cut -d/ -f 2 | cut -d. -f 1)
	
	singleton=$(echo $i | sed 's/R1//g')
	out_singleton=$(echo $nameR3 | cut -d/ -f 2 | cut -d. -f 1)
	
	samplename=$(echo $file | cut -d/ -f 2 | cut -d_ -f 1)
	
	## Trimmer les reads appariés à l’aide d’Alientrimmer
	java -jar soft/AlienTrimmer.jar -if $nameR1 -ir $nameR2 -c databases/contaminants.fasta -q 20 -of $2/trimmed/$out_nameR1.at.fq -or $2/trimmed/$out_nameR2.at.fq -os $2/trimmed/$out_singleton.at.sgl.fq
	
	## Reportez les distributions de qualités après filtrage qualité à l’aide de fastqc
	fastqc $2/trimmed/$out_nameR1.at.fq $2/trimmed/$out_nameR2.at.fq -o $2/fastqc_aprestrim
	
	## Fusionner les reads à l’aide Vsearch > Pairer les séquence F et R > Merger
	vsearch --fastq_mergepairs $2/trimmed/$out_nameR1.at.fq --reverse $2/trimmed/$out_nameR2.at.fq --fastaout $2/fasta/$samplename.fasta --label_suffix "$s;sample=$samplename;" --fastq_minovlen 40 --fastq_maxdiffs 15
done


#supprimer les espaces
#mettre tout dans les fastas dans un fichier amplicon
cat $2/fasta/*.fasta > $2/fasta/amplicon.fasta
sed -i -e "s/ //g" $2/fasta/amplicon.fasta


######## PARTIE 2 ########

#### Deduplication

	## Full length
soft/vsearch --derep_fulllength $2/fasta/amplicon.fasta --log=vsearch_fl_log --sizeout --minuniquesize 10 --output $2/fasta/amplicon_derep_fl.fasta --uc $2/fasta/amplicon_derep_fl.uc

	## Prefix
soft/vsearch --derep_prefix $2/fasta/amplicon.fasta --log=vsearch_pref.log --sizeout --minuniquesize 10 --output $2/fasta/amplicon_derep_pref.fasta --uc $2/fasta/amplicon_derep_pref.uc


#### Suppression des chimères
soft/vsearch --uchime_denovo $2/fasta/amplicon_derep_fl.fasta --nonchimeras $2/fasta/amplicon_non_chimeras.fasta --sizein --chimeras $2/fasta/chimeras.fasta


#### clustering
soft/vsearch --cluster_size $2/fasta/amplicon_non_chimeras.fasta --id 0.97 --sizein --sizeout --relabel OTU_ --centroids $2/fasta/OTU.fasta


######## PARTIE 3 ########

### Déterminer l'abondance des OTUs
soft/vsearch -usearch_global $2/fasta/amplicon.fasta --otutabout $2/fasta/OTU  -db $2/fasta/amplicon_dedfluchim_OTU.fasta --id 0.97


#### Annoter les OTU contre l’ensemble de 16S/18S fournis par le constructeur
soft/vsearch --usearch_global $2/fasta/OTU.fasta --userout  $2/fasta/annotations.fasta --db databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields query+target

sed '1iOTU\tAnnotation' -i $2/fasta/annotations.fasta

