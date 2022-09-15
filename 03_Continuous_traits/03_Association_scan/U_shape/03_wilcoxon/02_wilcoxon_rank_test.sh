#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=wilcoxon_signals
#SBATCH --nodes=1
#SBATCH --time=0-10:00:00
#SBATCH --partition normal
#SBATCH --cpus-per-task 1
#SBATCH --mem=50GB

#SBATCH --array=1-7
#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/log/02_wilcoxon_rank_test-%A_%a.out

###################
# Variables #
###################
echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

############ Run job ###########

file="/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/wilcoxon_test/input/list_traits.txt"

pheno=$(cat $file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/u_shape/06_wilcoxon/02_wilcoxon_rank_test.R $pheno

echo "Job done!"


