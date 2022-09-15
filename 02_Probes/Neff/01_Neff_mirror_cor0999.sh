#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=Neff_mirror   	    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-03:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=100GB           		# Memory required per node


#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/log/01_Neff_mirror_cor0999-%j.out

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Neff_after_22q112_definition/01_Neff_mirror_cor0999.R

echo "Job done"
