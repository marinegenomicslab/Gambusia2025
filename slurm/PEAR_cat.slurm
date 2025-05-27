#!/bin/bash

#SBATCH -J PEAR_cat                       # Name of the job
#SBATCH -o PEAR_cat.out                 # Name of file that will have program output
#SBATCH -e PEAR_cat.err                 # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                               # Partition
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=96:00:00

samtools merge -@64 -b bamlist.list -f cat-RRG.bam
samtools sort cat-RRG.bam -o cat-RRG.bam -@64
samtools index cat-RRG.bam
bedtools merge -i cat-RRG.bam -bed > mapped.bed

mawk -v OFS='\t' {'print $1,$2'} reference.fasta.fai > genome.file
cat namelist | parallel -j 48 "bedtools coverage -b {}-RG.bam -a mapped.bed -counts -sorted -g genome.file > {}.cov.stats"
