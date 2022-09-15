#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=prepare_wilcox_df     # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-01:00:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=50GB           			# Memory required per node


#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/log/01_prepare_wilcoxon_input_df-%j.out


#################
#    RUN JOB     #
##################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/u_shape/06_wilcoxon/01_prepare_wilcoxon_input_df.R

echo "Job Done!"
