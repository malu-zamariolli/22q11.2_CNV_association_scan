#####################################################
##### Combine results in the same file and delete individuals pheno files

#########################################
#### Packages

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom)

#########################################
#Combine results and exclude separate files 
list <- list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/", pattern = "_wilcoxon_results.txt", full.names = T, recursive = F)

print(paste0("There are ", length(list), " wilcoxon full results files. There should be 8 (1 per pheno)"))

df <- rbindlist(sapply(list, fread, simplify = F), use.names = T)
print(paste0("There are ", nrow(df), " rows in df"))

fwrite(df, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/8_traits_wilcoxon_all_results.txt", col.names = T, row.names = F, quote = F, sep = "\t")
unlink("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/*_wilcoxon_results.txt")

##################################################
#### Significant probes

list_sig <- list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/", pattern = "_wilcoxon_sig_results.txt", full.names = T, recursive = F)

print(paste0("There are ", length(list_sig), " wilcoxon sig results files"))

df_sig <- rbindlist(sapply(list_sig, fread, simplify = F), use.names = T)
print(paste0("There are ", nrow(df_sig), " rows in df"))

fwrite(df_sig, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/traits_sig_wilcoxon_results.txt", col.names = T, row.names = F, quote = F, sep = "\t")
unlink("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/*_wilcoxon_sig_results.txt")

#########################################
######## Add p-value to initial df (gwas results)

# Load CNV-GWAS dup-only results
gwas <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/gwas_results/continuous_trait22q112_sig_associations_dup.txt", header = T))

# Select columns of interest
df_wilcox <- df %>% select(phenotype, rsID, estimate, p.value)
colnames(df_wilcox)[c(3,4)] <- c("wilcoxon_estimate", "wilcoxon_pvalue")

# Merge with gwas results
merged <- merge(gwas, df_wilcox, by.y = c("phenotype", "rsID"), by.x = c("phenotype", "ID"))

if (nrow(merged) == nrow(gwas)) {
		print("All done!") 
	} else {
		print("Something is wrong :( ") 
}

# Save uptaded GWAS results
fwrite(merged, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/duplication_only/wilcoxon_test/dup_wilcoxon_GWAS_results.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


