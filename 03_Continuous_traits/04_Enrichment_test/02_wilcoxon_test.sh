#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=wilcoxon_neg_ctr   	    # Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-03:00:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=50GB           		# Memory required per node


#SBATCH --array=1-4
#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/log/02_wilcoxon_test_cont-%A_%a.out

###################
# Variables #
###################
echo "SLURM_JOB_ID: " $SLURM_JOB_ID
echo "SLURM_ARRAY_ID: " $SLURM_ARRAY_ID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID

############ Run job ###########

file="/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/path_list.txt"
input=$(awk '{print $1}' $file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
output=$(awk '{print $2}' $file | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

Rscript /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/11_traits_22q112_cnv_gwas/Continuous/Negative_controls/02_wilcoxon_test.R $input $output

echo "Job done!"


