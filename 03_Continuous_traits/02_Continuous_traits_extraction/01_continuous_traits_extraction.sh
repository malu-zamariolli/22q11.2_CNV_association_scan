#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=cont_traits  	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-01:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=10GB           		# Memory required per node

#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/log/03_continuous_traits_extraction-%j.out

#################
#   JOB INFO    #
#################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/09_traits_22q112/01_continuous_traits_extraction.R

echo "Job done"

