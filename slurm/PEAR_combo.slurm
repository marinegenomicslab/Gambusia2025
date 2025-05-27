#!/bin/bash

#SBATCH -J PEAR_combo                       # Name of the job
#SBATCH -o PEAR_combo.out                 # Name of file that will have program output
#SBATCH -e PEAR_combo.err                 # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                               # Partition
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=96:00:00

cat namelist | while read i; do
echo $i
echo -e "${i}_SE_filter.bam\n${i}_PE_filter.bam" > tmp.bamlist
samtools merge -h ${i}_PE_filter.bam -@64 -b tmp.bamlist -f $i-RG.bam
rm tmp.bamlist
samtools sort -@40 $i-RG.bam -o $i-RG.bam 2>> $i.log
samtools index ${i}-RG.bam
done
