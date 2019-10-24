#!/bin/bash
chmod -x 16S_workflow.sh

folderin=$1
folderout=$2

mkdir $2/fastqc_avanttrim
mkdir $2/fastqc_aprestrim
mkdir $2/trimmed
mkdir $2/vsearch

#### Deziper les fastq
gunzip $1/*.gz



#### Reportez les distributions de qualités avant et après filtrage qualité à l’aide de fastqc
fastqc $1/*_R1.fastq $1/*_R2.fastq -o $2/fastqc_avanttrim

for i in $list_fastq; do
    fastqc -o $2/fastqc_avanttrim --noextract fastq/$i
done

for file in $(ls $1/*_R1.fastq)
do
	nameR1="$file"
	nameR2=$(echo $file | sed "s/R1/R2/g")
	samplename=$(echo $file | cut -d/ -f 2 | cut -d_ -f 1)
	
	## Trimmer les reads appariés à l’aide d’Alientrimmer
	java -jar soft/AlienTrimmer.jar -if $nameR1 -ir $nameR2 -or $2/trimmed/$name_R1 -of $2/trimmed/$name_R2 -c databases/contaminants.fasta -q 20
	
	## fastqc
	fastqc $2/trimmed/*_R1.fastq.at.fq $2/trimmed/*_R2.fastq.at.fq -o $2/fastqc_aprestrim
	
	## Fusionner les reads à l’aide Vsearch > Pairer les séquence F et R > Merger
	vsearch --fastq_mergepairs $2/trimmed/*_R1.fastq.at.fq --reverse $2/trimmed/*_R2.fastq.at.fq --fastaout fasta/$samplename.fasta --label_suffix "$s;sample=$samplename;" --fastq_minovlen 40 --fastq_maxdiffs 15
done


#supprimer les espaces
#mettre tout dans les fastas dans un fichier amplicon
cat $2/vsearch/*.fasta > $2/vsearch/amplicon.fasta
sed -i -e "s/ //g" $2/vsearch/amplicon.fasta


#### Deduplication

	## Full length
soft/vsearch --derep_fulllength $2/vsearch/amplicon.fasta --log=vsearch_fl_log --sizeout --minuniquesize 10 --output $2/vsearch/amplicon_dedfl.fasta --uc $2/vsearch/amplicon_dedfl.uc

	## Prefix
soft/vsearch --derep_prefix $2/vsearch/amplicon.fasta --log=vsearch_pref.log --sizeout --minuniquesize 10 --output $2/vsearch/amplicon_dedpref.fasta --uc $2/vsearch/amplicon_dedpref.uc


#### Suppression des chimères
soft/vsearch --uchime_denovo $2/vsearch/amplicon_dedfl.fasta --nonchimeras $2/vsearch/amplicon_dedfluchim.fasta --sizein --chimeras $2/vsearch/chimeras_fl.fasta

soft/vsearch --uchime_denovo $2/vsearch/amplicon_dedpref.fasta --nonchimeras $2/vsearch/amplicon_dedprefuchim.fasta --sizein --chimeras $2/vsearch/chimeras_pref.fasta


#### clustering
soft/vsearch --cluster_size $2/vsearch/amplicon_dedfluchim.fasta --id 0.97 --sizein --sizeout --relabel OTU_ --centroids $2/vsearch/amplicon_dedfluchim_OTU.fasta

soft/vsearch --cluster_size $2/vsearch/amplicon_dedprefuchim.fasta --id 0.97 --sizein --sizeout --relabel OTU_ --centroids $2/vsearch/amplicon_dedprefuchim_OTU.fasta


### Déterminer l'abondance des OTUs
soft/vsearch -usearch_global $2/vsearch/amplicon.fasta --otutabout $2/vsearch/dedfluchim_OTU  -db $2/vsearch/amplicon_dedfluchim_OTU.fasta --id 0.97


soft/vsearch -usearch_global $2/vsearch/amplicon.fasta --otutabout $2/vsearch/dedprefuchim_OTU -db $2/vsearch/amplicon_dedprefuchim_OTU.fasta --id 0.97


#### Annoter les OTU contre l’ensemble de 16S/18S fournis par le constructeur
soft/vsearch --usearch_global $2/vsearch/amplicon_dedfluchim_OTU.fasta --userout  $2/vsearch/amplicon_dedfluchim_annot.fasta --db databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields query+target

soft/vsearch --usearch_global $2/vsearch/amplicon_dedprefuchim.fasta --userout  $2/vsearch/amplicon_dedprefuchim_annot.fasta --db databases/mock_16S_18S.fasta  --id 0.9 --top_hits_only --userfields query+target

sed '1iOTU\tAnnotation' -i $2/vsearch/amplicon_dedfluchim_annot.fasta
sed '1iOTU\tAnnotation' -i $2/vsearch/amplicon_dedprefuchim_annot.fasta

