# Calculate the CNV, duplication, and deletion frequency of all probes with genotype missingness < 5%

###############################################
### Libraries #################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#################################################
### STEP 1: Load genotype counts ################

# List all genotype count files from 01_genotype_count.sh
geno_count <- fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_count_chr22.frqx.gz", header = T)

print(paste0("Number of probes assessed on chr 22: ", nrow(geno_count)))


#################################################
### STEP 2: High-missingness probes #############
# Load high-missingness probes from 02_high_missingness_0.95.R
miss_probes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_missingness_0.95.txt", header = T, select = 2)) # loading rsIDs
print(paste0("Number of probes excluded due to high missingness: ", nrow(miss_probes)))

# Exclude high (5%) missingness probes
geno_count <- geno_count[!geno_count$SNP %in% miss_probes$SNP, ]
print(paste0("Number of probes remaining after excluding high missingness probes: ", nrow(geno_count)))

#####################################################################
### STEP 3: Convert genotype counts to number of CNV carriers #######

# Set up a dataframe
geno_freq <- data.frame(Chr = geno_count$CHR, SNP = geno_count$SNP, NumCNV = NA, NumDup = NA, NumDel = NA, NumCN = NA, NumMiss = NA)

# Number of duplication carriers (NumDup; > 2 copies) (Dup = TT)
geno_freq[which(geno_count$A1 == "A"), "NumDup"] <- geno_count[which(geno_count$A1 == "A"), "C(HOM A2)"] 
geno_freq[which(geno_count$A1 == "T"), "NumDup"] <- geno_count[which(geno_count$A1 == "T"), "C(HOM A1)"] 

# Number of deletion carriers (NumDel; < 2 copies) (Del = AA)
geno_freq[which(geno_count$A1 == "A"), "NumDel"] <- geno_count[which(geno_count$A1 == "A"), "C(HOM A1)"] 
geno_freq[which(geno_count$A1 == "T"), "NumDel"] <- geno_count[which(geno_count$A1 == "T"), "C(HOM A2)"]

# Number of CNV carriers (NumCNV; != 2 copies)
geno_freq$NumCNV <- geno_freq$NumDup + geno_freq$NumDel

# Number of copy-neutral carriers (NumCN, = 2 copies) (copy-neutral = AT)
geno_freq[, "NumCN"] <- geno_count[, "C(HET)"]

# Number of missing genotypes (NumMiss)
geno_freq[, "NumMiss"] <- geno_count[, "C(MISSING)"]


#################################################
### STEP 4: Genotype frequencies [%] ############

geno_freq$FreqCNV <- 100 * (geno_freq$NumCNV/(geno_freq$NumCNV + geno_freq$NumCN))
geno_freq$FreqDup <- 100 * (geno_freq$NumDup/(geno_freq$NumCNV + geno_freq$NumCN))
geno_freq$FreqDel <- 100 * (geno_freq$NumDel/(geno_freq$NumCNV + geno_freq$NumCN))


#################################################
### STEP 5: Print summary #######################

# CNVs
print("Number of probes with CNV frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqCNV >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqCNV >= 0.005), ])))

# Duplications
print("Number of probes with duplication frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqDup >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqDup >= 0.005), ])))

# Deletions
print("Number of probes with deletion frequency >= 0.005%:") 
print(as.data.frame(geno_freq[which(geno_freq$FreqDel >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq[which(geno_freq$FreqDel >= 0.005), ])))

######################################################
### STEP 6: Load probes for the 22q11.2 region and do summary

print("............Genotype (CNV) frequency for the 22q11.2 region (22:18400000-22500000).............")

# Probes come from .map file (/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Scripts_git/bepe_unil/03_CNV_GWAS_deletion_model/Neff/01_Probe_level_data.R)
probes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.map", header = F, select = c(2), col.names = c("rs"))) 
print(paste0("There are ", nrow(probes), " probes in the .map file for the 22q11.2 region"))

# Filter probes of interest:
geno_freq_region <- geno_freq[geno_freq$SNP %in% probes$rs, ]

# Summary
# CNVs
print("Number of probes with CNV frequency >= 0.005%:") 
print(as.data.frame(geno_freq_region[which(geno_freq_region$FreqCNV >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq_region[which(geno_freq_region$FreqCNV >= 0.005), ])))

# Duplications
print("Number of probes with duplication frequency >= 0.005%:") 
print(as.data.frame(geno_freq_region[which(geno_freq_region$FreqDup >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq_region[which(geno_freq_region$FreqDup >= 0.005), ])))

# Deletions
print("Number of probes with deletion frequency >= 0.005%:") 
print(as.data.frame(geno_freq_region[which(geno_freq_region$FreqDel >= 0.005), ] %>% count(Chr)))
print(paste0("TOTAL: ", nrow(geno_freq_region[which(geno_freq_region$FreqDel >= 0.005), ])))

#################################################
### Save ########################################

fwrite(geno_freq, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/genotype_frequency/genotype_freq_chr22_All.txt.gz", col.names = T, row.names = F, quote = F, sep = "\t" )

fwrite(geno_freq_region, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/genotype_frequency/genotype_freq_22:18400000_22500000_All.txt.gz", col.names = T, row.names = F, quote = F, sep = "\t" )

