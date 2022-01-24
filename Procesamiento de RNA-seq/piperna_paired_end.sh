#!/bin/bash
## This pipeline analysis RNA-seq data
## Author: Francisco J. Romero-Campero
## Contact: fran@us.es
if [ $# -eq 0 ]
  then
    echo "This pipeline analysis RNA-seq data."
    echo "Usage: piperna <param_file>"
    echo ""
    echo " param_file: File with the parameters specification. Please, check test/params.txt for an example."
    echo ""
    echo "Enjoy!"

    exit 0
fi


## Parameters loading
PARAMS=$1

WD=$(grep working_directory: $PARAMS | awk '{ print $2 }')
NUMSAM=$(grep number_of_samples: $PARAMS | awk '{ print $2 }')

SAMPLESLEFT=( )
SAMPLESRIGHT=( )

I=0
while [ $I -lt $NUMSAM ]
do
  SAMPLESLEFT[$I]=$(grep sample_$(($I + 1))_1: $PARAMS | awk '{ print $2 }')
  ((I++))
done

I=0
while [ $I -lt $NUMSAM ]
do
 SAMPLESRIGHT[$I]=$(grep sample_$(($I + 1))_2: $PARAMS | awk '{ print $2 }')
 ((I++))
done


## Debugging printing variable values
echo "Reading parameters from " $PARAMS 
echo "Input parameters: "
echo WD=$WD
echo NUMSAM=$NUMSAM

I=0
while [ $I -lt $NUMSAM ]
do
   echo sample_$((I+1))_1 = ${SAMPLESLEFT[$I]}
   ((I++))
done

I=0
while [ $I -lt $NUMSAM ]
do
   echo sample_$((I+1))_2 = ${SAMPLESRIGHT[$I]}
   ((I++))
done


## Generate working directory
mkdir $WD
cd $WD
mkdir genome annotation samples results logs
cd samples
I=1
while [ $I -le $NUMSAM ]
do
   mkdir sample$I
   ((I++))
done

## Generate genome index
cd $WD/genome
cp /mnt/space_158.42.124.151/home/frodmar/Mptak1_rna_processing/genome/* .

#cp /mnt/space_158.42.124.151/home/frodmar/data/genomes/Ostreococcus/GCF_000214015.3_version_140606_genomic.fa ./genome.fa
#gunzip genome.fa.gz

cd ../annotation
cp /mnt/space_158.42.124.151/home/frodmar/Mptak1_rna_processing/annotation/* .
#cp /mnt/space_158.42.124.151/home/frodmar/data/genomes/Ostreococcus/GCF_000214015.3_version_140606_genomic.gtf ./annotation.gtf
#wget -O annotation.gff https://marchantia.info/download/tak1v5.1/MpTak1v5.1_r1.gff
#gunzip annotation.gtf.gz

#cd ../genome

#hisat2-build genome.fa -p 6 index

## Copy samples
cd $WD/samples

I=0
while [ $I -lt $NUMSAM ] 
do
   #cp ${SAMPLES[$I]} sample$(($I+1))/sample$(($I+1)).fq.gz
   ln -s ${SAMPLESLEFT[$I]} sample$(($I+1))/sample$(($I+1))_1.fq.gz #JF211022
   ((I++))
done


I=0
while [ $I -lt $NUMSAM ]
do
   #cp ${SAMPLES[$I]} sample$(($I+1))/sample$(($I+1)).fq.gz
   ln -s ${SAMPLESRIGHT[$I]} sample$(($I+1))/sample$(($I+1))_2.fq.gz #JF211022
   ((I++))
done


I=1
while [ $I -le $NUMSAM ]
do 
   #sbatch -J sample$I -o $WD/logs/sample$I $HOME/opt/PIPERNA/rna_seq_sample_processing.sh $I $WD ${NUM_SAMPLES}
   bash $HOME/opt/piperna/rna_seq_sample_processing.sh $I $WD ${NUM_SAMPLES} #JF211022
   ((I++))
done
