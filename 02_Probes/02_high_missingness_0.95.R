# Probes with high genotype missingness (>5%) will be excluded 

###############################################
### Libraries #################################
library(data.table)
library(dplyr)


#################################################
### STEP 1: Load data ###################

# List all genotype count files from 01_genotype_count.sh
geno_count <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_count_chr22.frqx.gz"))

print(paste0("Number of probes assessed: ", nrow(geno_count)))

#################################################
### STEP 2: Detect high missingness probes ######

# Calculate the corresponding number of individuals (unrelated) for a missingness threshold at 5% 
Num <- nrow(as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_All.txt", header = F)))
missingness_thr <- Num*0.05

# Extract probes with high missingness
missingness_probes <- geno_count[which(geno_count$`C(MISSING)` > missingness_thr), ]
print(paste0("Number of probes with >5% missing genotypes in all individuals: ", nrow(missingness_probes)))
print(as.data.frame(missingness_probes %>% count(CHR)))


#################################################
### Save ########################################

fwrite(missingness_probes[,c(1,2)], "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/probes/genotype_missingness_0.95.txt", col.names = T, row.names = F, quote = F, sep = "\t")
