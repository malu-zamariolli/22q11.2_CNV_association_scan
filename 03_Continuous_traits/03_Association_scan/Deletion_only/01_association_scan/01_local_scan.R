# Local association study for the 22q11.2 region in the deletion-only model with traits linked to the 22q11.2 genes (HPO)

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

# Human Reference Build
print(paste0("Human reference build used: GRCh37/hg19"))

#################################################
### External arguments ##########################

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
start <- args[2]
end <- args[3]
print(paste0("Local association study for: chr", chr, ":", start, "-", end))


#################################################
### STEP 1: Temporary missingness file ##########

# Here we will perform a GWAS with all probes in the region of interest, only excluding high missingness probes (from script 02_high_missingness_0.95.R)

hm <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_missingness_0.95.txt"))
fwrite(data.frame(hm$SNP), "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/deletion_only/intermediate_files/high_missingness.txt", col.names = F, row.names = F, quote = F, sep = "\t")


#################################################
### STEP 2: Prepare PLINK input #################

# PLINK file set - deletion-only model
plink_file <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/deletion_only/ukb_cnv_chr", chr)

# INT covariate-corrected phenotypes (from 03_continuous_traits_extraction.R)
pheno_file <- "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_age_age2_sex_batch_PCs_All.txt"

# Unrelated samples (from 01_sample_filtering.R)
samples_file <- "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_All.txt"

# Probes to exclude due to high missingness (see STEP 1)
HM_file <- "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/deletion_only/intermediate_files/high_missingness.txt"

# Output file
output_file <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/deletion_only/intermediate_files/del", chr, ":", start, ":", end)

#################################################
### STEP 3: Perform local association study #####

# Run PLINK v2
system(paste("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/softwares/plink2/plink2",
 			 "--bfile", plink_file,
			 "--pheno",  pheno_file, "--no-psam-pheno",
			 "--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95", 
			 "--keep", samples_file,
			 "--exclude",  HM_file, "--chr", chr, "--from-bp", start, "--to-bp", end,
			 "--out", output_file))
unlink(HM_file)




