# Calculate the number of effective test, based on the method of Gao et al., 2008 (DOI: 10.1002/gepi.20310) 

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#############################################################################
### STEP 1: Load files

# data frame with gene coordinates
genes <- read.delim("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/Negative_controls/input/genes22q11.2_coordinates.txt")

# Genes selected for analysis
selected_genes <- fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/genes_selected_groups_continuous.txt")
selected_genes <- selected_genes$Gene
print(paste0(length(selected_genes), " genes will be analysed"))

## File with coordinate information for all probes
# All probes selected and cleaned in the 22q11.2 region; This corresponds to the .map file matching the binary PLINK file set for the mirror model (from 01_Probe_level_data.R)
probes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.map", header = F, select = c(1,2,4), col.names = c("chr", "rs", "pos")))

print(paste0("There are ", nrow(probes), " probes in the .map file for the 22q11.2 region"))

# File with probes used in analysis
probes_correlated <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Probes/864_probes_22q11.2.txt", header = T, select = c("ID"), col.names = "rs"))
print(paste0("All probes in the 22q11.2 region will be used: ", nrow(probes_correlated), " probes"))

#############################################################################
### STEP 2: Select probes per gene
# Obtain probes for gene of interest
# Get gene coordinates +- 10 kb
genes$Gene <- gsub(">", "", genes$Gene)
genes$Gene <- gsub("<", "", genes$Gene)
genes <- genes[genes$Gene %in% selected_genes, ]

# Add 10 kb intervals for each gene
genes$Start <- genes$Start - 10000
genes$End <- genes$End + 10000

# Select probes
probes_selected <- probes[probes$rs %in% probes_correlated$rs, ]
print(paste0(nrow(probes_selected), " will be used for gene selection"))

# Loop to obtain the probes inside the genes
datalist = list()

for (i in 1:nrow(genes)) {
  probes_temp <- probes_selected %>% filter(between(probes_selected$pos, 
                                           genes$Start[i], genes$End[i]))
  probes_temp$Gene <- genes$Gene[i]
  datalist[[i]] <- probes_temp 
}

probes_genes = do.call(rbind, datalist)

length(unique(probes_genes$rs)) 
rm(datalist)
rm(probes_temp)

#########################################################################
####### STEP 3: Calculate Number of Effective Tests for each gene

### 3A: Load the core files 

# All samples; This corresponds to the .ped file matching the binary PLINK file set 
all_samples <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.ped", header = F, select = 2, col.names = "eid"))
print(paste0("There are ", nrow(all_samples), " samples in the .ped file"))

# Samples to keep: are the results from sample_filtering.R
eids_to_keep <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_All.txt", header = F, col.names = "eid"))
print(paste0("There are ", nrow(eids_to_keep), " samples to keep"))


### 3B:  Loop for each gene's matrix #######

print("Starting loop to build each gene's matrix ...............")
genes_selected <- selected_genes

matrix_list <- list()
gene_probe <- list()
for (g in 1:length(genes_selected)) { 
	gene <- genes_selected[[g]]
	probes_to_keep <- probes_genes[probes_genes$Gene %in% gene, "rs"]
  # Save number of probes per gene
	df <- data.frame("Gene" =  gene, "Number_of_probes:" = length(unique(probes_to_keep)))
	gene_probe[[g]] <- df
  # Load .ped chunk
	cnv <- matrix()
	ped_chunck <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/plink_input/probe_level_mirror/ukb_cnv_22:18400000_22500000.ped", header = F, drop = c(1:6)))
  
	rownames(ped_chunck) <- all_samples$eid
	ped_chunck <- ped_chunck[rownames(ped_chunck) %in% eids_to_keep$eid, ]
  
	cnv_chunck <- data.frame(mapply(paste0, ped_chunck[][c(T, F)], ped_chunck[][c(F,T)]), stringsAsFactors = F)
	rownames(cnv_chunck) <- rownames(ped_chunck); rm (ped_chunck)
	colnames(cnv_chunck) <- probes[ , "rs"] 
	cnv_chunck <- cnv_chunck[, colnames(cnv_chunck) %in% probes_to_keep]
 	cnv <- cnv_chunck
	matrix_list[[g]] <- cnv
	names(matrix_list)[g] <- gene
	print(paste0("Final dimensions for the ", gene, " : ", nrow(cnv), " (samples; row) x ", ncol(cnv), " (probes; column)"))
}

all_gene_probe = do.call(rbind, gene_probe)
fwrite(all_gene_probe, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/selected_genes_number_probes.txt", sep = "\t", row.names = F, col.names = T, quote = F)

### 3C: Calculate Neff

print("Starting loop to calculate each gene's Neff ..............")

Neff_list <- list()
for (i in 1:length(matrix_list)) {
  matrix_list[[i]][matrix_list[[i]] == "AA"] <- -1
  matrix_list[[i]][matrix_list[[i]] == "AT" | matrix_list[[i]] == "TA"] <- 0
  matrix_list[[i]][matrix_list[[i]] == "TT"] <- 1
  matrix_list[[i]][matrix_list[[i]] == "00"] <- NA
  matrix_list[[i]] <- as.data.frame(matrix_list[[i]])
  matrix_list[[i]] <- mutate_all(matrix_list[[i]], function(x) as.numeric(as.character(x)))
  # Tranform in matrix
  matrix_list[[i]] <- as.matrix(matrix_list[[i]])
  # Calculate the correlation matrix
  cnv_cor <- cor(matrix_list[[i]], use = 'pairwise.complete.obs')
  # Set missing values to 0 to allow svd()
  cnv_cor[which(is.na(cnv_cor))] <- 0
  #### Calculate the eigenvalues
  # svd() computes the singular-value decomposition of a recatngular matrix, with $d being a vector containing the singular values of the decomposition
  cnv_EV <- svd(cnv_cor)$d 
  #### Calculate Neff - # Neff is defined as the number of eigenvalues required to explain 99.5% of the variation 
	sum_EV <- 0
	count <- 0
	while(sum_EV/sum(cnv_EV) < 0.995) {
		count <- count + 1
		sum_EV <- sum_EV + cnv_EV[count]
	}
  Neff_list[[i]] <- data.frame("Neff" = count, "Gene" = c(names(matrix_list)[i]))
  print(paste0("Neff : ", count))
}

all_Neff = do.call(rbind, Neff_list)

# Merge with number of probes
merged <- merge(all_Neff, all_gene_probe, by = "Gene")

fwrite(merged, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Negative_controls/continuous/Neff_gene.txt", col.names = T, row.names = F, quote = F, sep = "\t")




