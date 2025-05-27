#!/bin/bash

#SBATCH -J PEAR_PE_filter                     # Name of the job
#SBATCH -o PEAR_PE_filter.out                 # Name of file that will have program output
#SBATCH -e PEAR_PE_filter.err                 # Name of the file that will have job errors, if any
#SBATCH -N 1                                  # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                                 # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                             # Partition
#SBATCH --mail-type=begin                     # email me when the job starts
#SBATCH --mail-type=end                       # email me when the job ends
#SBATCH --time=96:00:00

cat namelist | xargs -I{} -P48 $WORK/slurm/pear/bam_PE_filter.sh {}
