#!/bin/bash
#SBATCH --job-name=rna_seq_sample_processing_fernando
#SBATCH --export=ALL
#SBATCH --output=rna_seq_sample_processing

## Author: Francisco J. Romero-Campero
## Contact: fran@us.es

## Reading input parameters
SAMPLE_ID=$1
WD=$2
NUM_SAMPLES=$3

## Access sample folder
cd $WD/samples/sample${SAMPLE_ID}

## QC
fastqc sample${SAMPLE_ID}_1.fq.gz -t 6
fastqc sample${SAMPLE_ID}_2.fq.gz -t 6

## Mapping to reference genome
hisat2 --dta -x $WD/genome/index -1 sample${SAMPLE_ID}_1.fq.gz -2 sample${SAMPLE_ID}_2.fq.gz -p 6 -S sample${SAMPLE_ID}.sam 
samtools sort -o sample${SAMPLE_ID}.bam --threads 6 sample${SAMPLE_ID}.sam

## Transcript assembly
stringtie -G $WD/annotation/annotation.gff -p 6 -o sample${SAMPLE_ID}.gtf -l sample${SAMPLE_ID} sample${SAMPLE_ID}.bam
rm sample${SAMPLE_ID}.sam

##Quantification
stringtie -e -B -G $WD/annotation/annotation.gff -p 6 -o sample${SAMPLE_ID}.gtf -A gene_abun${SAMPLE_ID}.tab sample${SAMPLE_ID}.bam

## Write information about sample transcriptome
echo $WD/samples/sample${SAMPLE_ID}/sample${SAMPLE_ID}.gtf >> $WD/results/merge_list.txt

## Synchronization point through blackboards
echo "sample${SAMPLE_ID} DONE" >> $WD/logs/blackboard 

DONE_SAMPLES=$(wc -l $WD/logs/blackboard | awk '{ print $1 }')

#if [ ${DONE_SAMPLES} -eq ${NUM_SAMPLES} ]
#then
    #sbatch -J transcriptome_merging -o $WD/logs/transcriptome_merging $HOME/opt/PIPERNA/transcriptome_merging.sh $WD 
#    $HOME/opt/PIPERNA/transcriptome_merging.sh $WD #JF211022
#fi
