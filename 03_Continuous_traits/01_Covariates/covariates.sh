#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=cov  	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-01:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=10GB           		# Memory required per node

#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/log/covariates-%j.out

#################
#   JOB INFO    #
#################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Scripts/covariates.R

echo "Job done"

