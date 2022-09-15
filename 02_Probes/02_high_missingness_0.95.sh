#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=geno_miss    	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:10:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           		# Memory required per node

#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/log/02_high_missingness-%j.out


#################
#   JOB INFO    #
#################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/scripts/01_frequency/02_high_missingness_0.95.R 

echo "Job done"

