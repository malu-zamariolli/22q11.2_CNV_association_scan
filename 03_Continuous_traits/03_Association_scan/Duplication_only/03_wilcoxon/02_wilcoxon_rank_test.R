# Perform wilcoxon rank test for significant continuous traits to verify CNV-GWAS results

#######################################
####### Packages 

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom)


#######################################
# External arguments
cmd_args = commandArgs(trailingOnly=TRUE)
pheno = cmd_args[1]

print(pheno)

########################################
#### Load data

# Load CNV-GWAS dup-only model results
gwas <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/gwas_results/continuous_trait22q112_sig_associations_dup.txt", header = T))

# Vector with signals (of interest) (column - ID)
probes <- as.vector(unique(gwas$ID))

# Load input df for wilcoxon (script 01_prepare_wilcoxon_input_df.R)
df <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/intermediate_files/significant_continuousTraits_probeLevel_raw.txt", header = T))


#########################################
####### Wilcoxon rank test
print(paste0("Starting wilcoxon rank test for pheno: ........", pheno))

# Wilcoxon rank test : wilcox.test(continuous ~ binary, data, conf.int = T)

df[, probes] <- lapply(df[, probes], as.factor)

snp <- unique(gwas[gwas$phenotype == pheno, "ID"])
print(paste0("Number of SNPs to be tested for ",  pheno, ":" ,length(snp)))

wilcox <- map(df[,snp], ~wilcox.test(df[,pheno] ~ .x, conf.int = T, na.action = na.omit)) %>% map_dfr(broom::tidy, .id = "rsID")
wilcox$phenotype <- pheno


fwrite(wilcox, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/", pheno, "_dup_wilcoxon_results.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

print(paste0("Method: Wilcoxon rank sum test with continuity correction. Alternative: two.sided"))

thr <- 0.01

wilcox_sig <- wilcox %>% filter(p.value < thr)
print(paste0("There are ", nrow(wilcox_sig), " significant results from wilcoxon test for pheno: ", pheno))
print(paste0("There were ", length(snp), " significant results from CNV-GWAS"))

fwrite(wilcox_sig, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/", pheno, "_dup_wilcoxon_sig_results.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


