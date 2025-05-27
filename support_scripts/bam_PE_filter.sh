#!/bin/bash

#Function to filter PE bam files
#Usage: bam_PE_filter.sh <sample_name>

#Getting Sample name in the variable list
i=$1

#Noting progress
#echo "Filtering $i"

#Checking for PE file
if [ ! -s ${i}_PE.bam ]; then exit; fi

#Filter for PE reads on the same contig or SE reads
#Remove reads with mate unmapped
samtools view -hq 40 ${i}_PE.bam | awk '$7 == "=" || $1 ~ /^@/' | samtools view -bS - | samtools view -hF 0x8 - |  samtools view -bS - > ${i}_tmp.bam       #https://www.biostars.org/p/118301/

#Finding reads flagged as not properly paired
samtools view -F 0x2 ${i}_tmp.bam > ${i}_exclude.list

#Removing not properly paired if appropriate
TMP=$(cat ${i}_exclude.list | wc -l)
TOT=$(samtools view -c ${i}_tmp.bam)
if [ $TOT -eq 0 ]; then PER=0
else PER=$(echo "($TMP*1000)/$TOT" | bc)
fi

echo "$TMP reads out of $TOT found to be not properly paired in $i"

if [ $TMP -eq 0 ]; then mv ${i}_tmp.bam ${i}_PE_filter.bam
rm ${i}_exclude.list
elif [ $PER -lt 10 ]; then
/work/marinegenomics/afields3/programs/miniconda3/bin/java -jar /work/marinegenomics/afields3/programs/miniconda3/share/picard-3.1.1-0/picard.jar FilterSamReads I=${i}_tmp.bam O=${i}_PE_filter.bam READ_LIST_FILE=${i}_exclude.list FILTER=excludeReadList  2> ${i}_PE_fil.log
rm ${i}_tmp.bam ${i}_exclude.list
else echo $i >> Check_Samples.txt
/work/marinegenomics/afields3/programs/miniconda3/bin/java -jar /work/marinegenomics/afields3/programs/miniconda3/share/picard-3.1.1-0/picard.jar FilterSamReads I=${i}_tmp.bam O=${i}_PE_filter.bam READ_LIST_FILE=${i}_exclude.list FILTER=excludeReadList  2> ${i}_PE_fil.log
rm ${i}_tmp.bam ${i}_exclude.list
fi
