#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=gwas_merge_del       # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:40:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=1GB           			# Memory required per node

#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/deletion_only/log/03_merge_del-%j.out


##################
#    RUN JOB     #
##################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/deletion_only/03_gwas/03_merge_del.R

echo "Job Done!"
