#!/bin/bash

#SBATCH -J PEAR_trimmed                       # Name of the job
#SBATCH -o PEAR_trimmed.out                   # Name of file that will have program output
#SBATCH -e PEAR_trimmed.err                   # Name of the file that will have job errors, if any
#SBATCH -N 1                                  # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                             # Partition
#SBATCH --mail-type=begin                     # email me when the job starts
#SBATCH --mail-type=end                       # email me when the job ends
#SBATCH --time=96:00:00

ls *.R1.fq.gz | sed 's/.R1.fq.gz//g' > namelist
mkdir logfiles

cat namelist | while read i; do
echo "Processing $i"
pearRM -f $i.R1.fq.gz -r $i.R2.fq.gz -m 290 -n 97 -o $i -j 64 1> logfiles/${i}_pear.log
done
