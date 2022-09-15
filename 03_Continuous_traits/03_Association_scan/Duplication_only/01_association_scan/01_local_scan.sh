#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=cnv_gwas_cont_dup			# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          				# Define number of nodes required
#SBATCH --time=0-00:40:00   			# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  			# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=3GB           			# Memory required per node

#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/log/01_local_scan-%j.out


##################
#   VARIABLES    #
##################

chr=22                                                  # Chromosome
start=18400000                                         # Start position
end=22500000                                            # End position



##################
#    RUN JOB     #
##################

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/duplication_only/03_gwas/01_local_scan.R ${chr} ${start} ${end} ${pheno}

echo "Job Done!"
