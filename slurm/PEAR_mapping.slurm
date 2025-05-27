#!/bin/bash

#SBATCH -J PEAR_mapping                      # Name of the job
#SBATCH -o PEAR_mapping.out                  # Name of file that will have program output
#SBATCH -e PEAR_mapping.err                  # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                     # Partition
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=96:00:00

cat namelist | while read i; do
echo $i
bwa mem -L 20,5 -t 60 -a -M -T 20 -A 1 -B 3 -O 5 -R "@RG\tID:$i\tSM:$i\tPL:Illumina" reference.fasta $i.unassembled.forward.fastq $i.unassembled.reverse.fastq  2> bwa.${i}_PE.log | mawk '$6 !~/[2-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -@60 -q 1 -SbT reference.fasta - > ${i}_PE.bam 2>${i}_PE.bam.log &
bwa mem -L 20,5 -t 60 -a -M -T 20 -A 1 -B 3 -O 5 -R "@RG\tID:$i\tSM:$i\tPL:Illumina" reference.fasta $i.assembled.fastq 2> bwa.${i}_SE.log | mawk '$6 !~/[2-9].[SH]/ && $6 !~ /[1-9][0-9].[SH]/' | samtools view -@60 -q 1 -SbT reference.fasta - > ${i}_SE.bam 2>${i}_SE.bam.log
samtools view -hBq 40 -@ 64 ${i}_SE.bam > ${i}_SE_filter.bam
wait
done
