#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=gwas_correct_dup     # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-01:00:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=3GB           			# Memory required per node


#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/log/02_correct_dup-%j.out


#################
#    RUN JOB     #
##################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/duplication_only/03_gwas/02_correct_dup.R

echo "Job Done!"
