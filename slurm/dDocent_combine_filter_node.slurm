#!/bin/bash

#SBATCH -J dDocent_combine_f                      # Name of the job
#SBATCH -o dDocent_combine_f.out                  # Name of file that will have program output
#SBATCH -e dDocent_combine_F.err                  # Name of the file that will have job errors, if any
#SBATCH -N 1                                 # Number of nodes ( the normal cluster partion has 22 total )
#SBATCH -n 5                                # Number of cores ( my test allocated 2 per node )
#SBATCH -p jgoldq,normal                               # Partition
#SBATCH --mail-type=begin                    # email me when the job starts
#SBATCH --mail-type=end                      # email me when the job ends
#SBATCH --time=24:00:00

#module load ddocent/2.9.1
vcfcombine SNP.raw_*.QC.recode.vcf | sed -e 's/       \.\:/   \.\/\.\:/g' > cat_f.vcf
