### Prepare dataframe that will be input for wilcoxon test

####### Packages 

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom)

#######################

#### Load data

# Load CNV-GWAS dup-only model results (full table with annotated frequencies)
gwas <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/gwas_results/continuous_trait22q112_sig_associations_dup.txt", header = T))

# Vector with signals (of interest) (column - ID)
probes <- as.vector(unique(gwas$ID))

# Vector with traits (column PHENO)
traits_vector <- as.vector(unique(gwas$phenotype))

##############################################################
######### Prepare data frame with probe_level data for duplication x CN (deletion will be considered NA)

# Load probe-level data 
probe_level <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Association_22region_binary/probe_level_input/probe_level_cn_22:18400000_22500000.txt", header = T))
print(paste0("There are ", nrow(probe_level), " rows in probe level df"))
print(paste0("There are ", ncol(probe_level), " cols in probe level df"))

# AA -> -1 (= deletion); AT -> 0 (= copy-neutral); TT -> 1 (= duplication)
# Deletion will be assigned as missing values (NA)
probe_level[,-1][probe_level[,-1] == -1] <- NA


fwrite(probe_level, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/probe_level_data/22q11_2_copynumber_status_dup_only_model.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

#########################################################################

# Probe-level data (only eid + probe columns of interest)
probe_data <- probe_level[, c("eid", probes)]
print(paste0("There are ", nrow(probe_data), " individuals in the probe level data df"))
print(paste0("There are ", ncol(probe_data), " probes (SNPs) in the probe level data df"))

# Load phenotype data (continuous, corrected) only for trait of interest (significant traits)
traits <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_age_age2_sex_batch_PCs_All.txt", header = T, select = c("IID", traits_vector)))
print(paste0("There are ", nrow(traits), " individuals in the traits data df"))
print(paste0("There are ", ncol(traits) -1, " traits in the traits data df"))
colnames(traits)[[1]] <- "eid"

##########################################
####### Merge probe-level and traits 
# Merge probe-level and traits data by eid

df <- merge(traits, probe_data, by = "eid")
print(paste0("There are ", nrow(df), " individuals in merged df. There should be: ", nrow(probe_data)))

	## Save as intermediate file
fwrite(df, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/intermediate_files/significant_continuousTraits_probeLevel_raw.txt", col.names = T, row.names = F, quote = F, sep = "\t")

#########################################

