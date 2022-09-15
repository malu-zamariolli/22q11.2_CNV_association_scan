#!/bin/bash

#################
#   RUN INFO    #
#################

#SBATCH --job-name=geno_count   	# Give your job a name to recognize it in queue overview
#SBATCH --nodes=1          			# Define number of nodes required
#SBATCH --time=0-00:40:00   		# Define how long the job will run (max. is 7 days, default is 1h)
#SBATCH --partition=normal  		# Define the partition on which teh job shall run. May be omitted
#SBATCH --mem=3GB           		# Memory required per node


#SBATCH --output=/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/log/01_genotype_count-%j.out


# Prepare input 
input=$(echo "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/mirror/ukb_cnv_chr22") # This corresponds to the .bim of the PLINK file set described in the methods of the paper (Table 1)
samples_all=$(echo "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_plink1.9_All.txt") # results from /01_samples/sample_filtering.R
output_all=$(echo "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_count_chr22")


#################
#   JOB INFO    #
#################

# PLINK v1.9 
/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/softwares/plink/plink \
--bfile ${input} \
--keep ${samples_all} \
--freqx gz \
--out ${output_all}


echo "Job done"


