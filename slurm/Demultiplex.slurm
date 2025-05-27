#!/bin/bash

#SBATCH -J demultiplex          # Name of the job
#SBATCH -o demultiplex.out      # Name of file that will have program output
#SBATCH -e demultiplex.err      # Name of the file that will have job errors, if any
#SBATCH -N 1                    # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 64                   # Number of cores ( my test allocated 2 per node )
#SBATCH -p normal               # Partition
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         #email me when the job ends
#SBATCH --time=24:00:00

module load stacks/2.66
sh $FILE
