# Calculate the number of effective test for the 22q11.2 region, based on the method of Gao et al., 2008 (DOI: 10.1002/gepi.20310) 
# The number of effective tests for the mirror model will be more stringent than for the deletion-only or duplication-only models separately

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Clean the .ped file #################

### Load the core files #########################

# All samples; This corresponds to the .ped file matching the binary PLINK file set for the mirror model 

all_samples <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.ped", header = F, select = 2, col.names = "eid"))
print(paste0("There are ", nrow(all_samples), " samples in the .ped file"))

# All probes selected and cleaned in the 22q11.2 region; This corresponds to the .map file matching the binary PLINK file set for the mirror model 
all_probes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.map", header = F, select = c(1,2), col.names = c("chr","rs"))) 
print(paste0("There are ", nrow(all_probes), " probes in the .map file for the 22q11.2 region"))

# load probes to keep (from correlation)
probes_to_keep <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Probes/504_probes_r2_999.txt", header = F, col.names = "rs"))
print(paste0("There are ", nrow(probes_to_keep), " probes to keep (correlated at r2 0.999 with at least 10 other probes (LD-friends)) "))


### Load .ped  & clean it #######

cnv <- matrix()

ped_chunck <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.ped", header = F, drop = c(1:6)))
	
# Clean samples by retaining unrelated samples (rows)
rownames(ped_chunck) <- all_samples$eid

# Clean variants by retaining filtered and pruned probes (columns)
cnv_chunck <- data.frame(mapply(paste0, ped_chunck[][c(T, F)], ped_chunck[][c(F,T)]), stringsAsFactors = F)


rownames(cnv_chunck) <- rownames(ped_chunck); rm (ped_chunck)
colnames(cnv_chunck) <- all_probes[ , "rs"]



cnv <- cnv_chunck[, colnames(cnv_chunck) %in% probes_to_keep$rs]
	
print(paste0("Final dimensions for the cleaned cnv file: ", nrow(cnv), " (samples; row) x ", ncol(cnv), " (probes; column) for r^2 "))

### Recode & make numeric matrix ################

# AA -> -1 (= deletion); AT -> 0 (= copy-neutral); TT -> 1 (= duplication)
cnv[cnv == "AA"] <- -1
cnv[cnv == "AT" | cnv == "TA"] <- 0
cnv[cnv == "TT"] <- 1
cnv[cnv == "00"] <- NA
cnv <- mutate_all(cnv, function(x) as.numeric(as.character(x)))

# Tranform in matrix
cnv <- as.matrix(cnv)

print(paste0("Dimensions of matrix :", dim(cnv)))

# STEP 2: Calculate the correlation matrix ###

# Calculate the correlation matrix
cnv_cor <- cor(cnv, use = 'pairwise.complete.obs')
rm(cnv)

# Set missing values to 0 to allow svd()
cnv_cor[which(is.na(cnv_cor))] <- 0

# STEP 3: Calculate the eigenvalues ###########
# svd() computes the singular-value decomposition of a rectangular matrix, with $d being a vector containing the singular values of the decomposition
cnv_EV <- svd(cnv_cor)$d 

# STEP 4: Calculate Neff ######################

# Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the CNV data 
sum_EV <- 0
count <- 0

while(sum_EV/sum(cnv_EV) < 0.995) {
	count <- count + 1
	sum_EV <- sum_EV + cnv_EV[count]
}

print(paste0("Neff on 22q11.2 region for probes correlated at r^2 0.999 with at least 10 other probes : ", count))
# Save
fwrite(data.frame(count), paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Neff/Neff_r2_999_22_18630000_21910000.txt"), col.names = F, row.names = F, quote = F, sep = "\t")






