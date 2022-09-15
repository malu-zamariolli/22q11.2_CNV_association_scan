############ Wilcoxon rank test - one sided 
## Compare associations between group of linked traits and non-linked traits (negative controls) per gene

# Wilcoxon rank test : wilcox.test(continuous ~ binary, data, conf.int = T)
  # binary : group (associated or negative)
  # continuous - p-value

###############
### Packages
library(tidyverse)
library(data.table)
#####################

###########################################################################################
####### STEP 1: Load files ##############################################################

cmd_args = commandArgs(trailingOnly=TRUE)
input = cmd_args[1]
output = cmd_args[2]

print(input)
 
print("Loading files ....")

## File with genes and Neff from each
gene_neff <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/Neff_gene.txt", header = T))

## File with trait-gene
gene_trait <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/18continuous_traits_HPO_gene.txt", header = T))
gene_trait <- gene_trait %>% select(Pheno, Gene)
gene_trait <- unique(gene_trait)

## Filter for genes that have been selected based on number of traits per group 
selected_genes <- as.vector(gene_neff$Gene)
print(paste0(length(selected_genes), " genes are being tested"))

### data frame with gene coordinates
genes <- read.delim("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/Negative_controls/input/genes22q11.2_coordinates.txt")

## File with coordinate information for all probes
# All probes selected and cleaned in the 22q11.2 region; This corresponds to the .map file matching the binary PLINK file set for the mirror model (from 01_Probe_level_data.R)
probes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.map", header = F, select = c(1,2,4), col.names = c("chr", "rs", "pos")))

# Load probes of interest
probes_to_keep <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Probes/864_probes_22q11.2.txt", header = T,select = c("ID"), col.names = "rs"))
print(paste0("For wilcoxon, all probes in the 22q11.2 region will be used: ", nrow(probes_to_keep), " probes"))

### Load association results for each model
df <- as.data.frame(fread(paste0(input), header = T))
df <- df[df$ID %in% probes_to_keep$rs, ]
	 
df <- df %>% select(phenotype, ID, P)

###################################################################################
############# STEP 2: Obtain probes for genes of interest ########################
print("Obtain probes per gene .....")

probes <- probes[probes$rs %in% probes_to_keep$rs, ]
print(paste0(nrow(probes), " will be used for gene selection"))

# Get gene coordinates +- 10 kb
genes$Gene <- gsub(">", "", genes$Gene)
genes$Gene <- gsub("<", "", genes$Gene)
genes <- genes[genes$Gene %in% selected_genes, ]

# Add 10 kb intervals from each gene
genes$Start <- genes$Start - 10000
genes$End <- genes$End + 10000

# Loop to obtain the probes inside the genes
datalist = list()

for (i in 1:nrow(genes)) {
  probes_temp <- probes %>% filter(between(probes$pos, 
                                           genes$Start[i], genes$End[i]))
  probes_temp$Gene <- genes$Gene[i]
  
  datalist[[i]] <- probes_temp # add it to your list
}

probes_genes = do.call(rbind, datalist)

length(unique(probes_genes$rs)) 
rm(datalist)
rm(probes_temp)

###################################################################################
#### STEP 3: Create input df for wilcoxon
## Data frame will have traits in rows, genes with group info in columns and a column per probe with p-values

print("Creating input for wilcoxon test .......")

# Create matrix (df) with traits in rows and Genes in columns
mat_gene_trait <- as.data.frame.matrix(table(gene_trait))
mat_gene_trait <- mat_gene_trait[, c(selected_genes)]
		# 1 will be associated traits group and 0 will be negative control group

# Pivot association results and do df for wilcoxon test
df <- df %>% pivot_wider(names_from = ID, values_from = P)

## Combine p-values df with Gene groups
mat_gene_trait <- cbind(phenotype = rownames(mat_gene_trait), data.frame(mat_gene_trait, row.names=NULL))
mat_gene_trait[, c(-1)] <- lapply(mat_gene_trait[, c(-1)], as.factor)

## Merge dfs
df_group <- merge(mat_gene_trait, df, by = "phenotype")

#####################################################################################
###########	STEP 4: Wilcoxon test #############################################

print("Starting wilcoxon test  .........")

list_wilcox = list()
sig_wilcox = list()
for (i in 1:length(selected_genes)) {
  i2 <- selected_genes[i]
  probes_i <- probes_genes %>% filter(Gene == i2)
  probes_i <- probes_i$rs
  thr <- gene_neff[gene_neff$Gene == i2, "Neff"]
  datalist = list()
  for (p in 1:length(probes_i)) {
    p2 <- probes_i[p]
	print(paste0("Testing gene: ", i2, "and probe: ", p2)) # Added this to identify possible errors
	# Do wilcoxon with tryCatch to bypass error (for some probes there are no association values for both groups)
	wilcox_temp <- tryCatch(broom::tidy(wilcox.test(df_group[ ,p2] ~ df_group[ ,i2], data = df_group, conf.int = T, na.action = na.omit, alternative = "greater")), error = function(e) return(data.frame("estimate" = NA, "statistic" = NA, "p.value" = NA, "conf.low" = NA, "conf.high" = NA, "method" = NA, "alternative" = NA,"Gene" = NA, "probe" = NA)))  
    wilcox_temp$Gene <- i2
    wilcox_temp$probe <- p2
    datalist[[p]] <- wilcox_temp
  }
  wilcox = do.call(rbind, datalist)
  rm(datalist)
  rm(wilcox_temp)
  list_wilcox[[i]] <- wilcox
  print("Selecting significant signals")
  wilcox_sig <- wilcox %>% filter(wilcox$p.value < 0.05/thr)
  sig_wilcox[[i]] <- wilcox_sig

}

wilcox_all = do.call(rbind, list_wilcox)
fwrite(wilcox_all, paste0(output, "wilcoxon_all_per_probe_continuous.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
wilcox_sig = do.call(rbind, sig_wilcox)
fwrite(wilcox_sig, paste0(output, "wilcoxon_sig_per_probe_continuous.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


rm(i)
print(paste0(length(unique(wilcox_sig$Gene)), " genes are significant out of ", length(selected_genes), " genes tested")) 

### Check if the smallest p-value of each gene is significant
print("Selecting smallest significant p-values ..............")

df <- wilcox_all %>% group_by(Gene) %>% slice(which.min(p.value))
df <- df %>% select(Gene, probe, p.value)
df$thr <- NA

# Vector with genes that had wilcox results
genes <- as.vector(df$Gene)

list_df <- list()
for (i in 1:length(genes)) {
  i2 <- genes[i]
  thr <- gene_neff[gene_neff$Gene == i2, "Neff"]
  df[df$Gene == i2, "thr"] <- 0.05/thr
  p <- df[df$Gene == i2, "p.value"]
  t <- df[df$Gene == i2, "thr"]
  if (p < t) {
  	list_df[[i]] <- df[df$Gene == i2, ]
     } else { print(paste0("Gene not significant: ", i2))
  }
}

df_all <- do.call(rbind, list_df)

# Save data frame with probes with smallest p-value 
fwrite(df, paste0(output, "wilcoxon_smallest_p_continuous.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

if (length(list_df) > 0) {
	fwrite(df_all, paste0(output, "wilcoxon_sig_smallest_p_continuous.txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
	
		} else { print("No significant genes")
}


rm(list_wilcox)
rm(list_df)

