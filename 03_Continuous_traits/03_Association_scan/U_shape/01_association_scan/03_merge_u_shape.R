# Merge files by phenotype and extract significant associations 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(gtools)

# ************************************************ GWAS RESULTS ************************************************ #
print("Merging GWAS data")
# Define threshold for significance (FWER 5%). 
thr <- 0.05/(6*(16+113))
print(paste0("p-value: ", thr))

##########################################
### STEP 1: List phenotypes ##############

# These are the corrected output files from the GWAS
list_files <- list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/intermediate_files/", pattern = ".glm.linear", full.names = T, recursive = F)
pheno_list <- unique(gsub(".*\\.", "", sub(".glm.linear", "", list_files)))
print(paste0("There are ", length(pheno_list), " unique phenotypes."))


#########################################
### STEP 2: Merge phenotypes in the same file
all_pheno <- rbindlist(sapply(list_files, fread, simplify = F), use.names = T, idcol = 'file')
all_pheno$file <- gsub("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/intermediate_files//u_shape22:18400000:22500000.", "", all_pheno$file)
all_pheno$file <- gsub(".glm.linear", "", all_pheno$file)

all_pheno <- all_pheno[order(-P),]

colnames(all_pheno)[1] <- "phenotype"

fwrite(all_pheno, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/gwas_results/continuous_trait22q112_associations_u_shape.txt", col.names = T, row.names = F, quote = F, sep = "\t")

##########################################
### STEP 3: Significant hits ############

sig_all <- all_pheno[which(all_pheno$P <= thr), ]
print(paste0(nrow(sig_all), " significant associations for phenotypes"))
print("Significant phenotypes: ")
print(unique(sig_all$phenotype))
print(paste0(print(length((unique(sig_all$phenotype)))), " out of 18 phenotypes are significant"))

fwrite(sig_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/gwas_results/continuous_trait22q112_sig_associations_u_shape.txt", col.names = T, row.names = F, quote = F, sep = "\t")

################################
### STEP 4: p < 0.05 hits

nominal <- all_pheno[which(all_pheno$P <= 0.05), ]
print(paste0(nrow(nominal), " nominal significance"))
print("Nominal significance phenotypes: ")
print(unique(nominal$phenotype))
print(paste0(print(length((unique(nominal$phenotype)))), " out of 18 phenotypes"))

fwrite(nominal, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/gwas_results/continuous_trait22q112_nominal_associations_u_shape.txt", col.names = T, row.names = F, quote = F, sep = "\t")


