# Calculate the number of effective test for the traits being evaluated

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#################################################
### STEP 1: Load files #################


# Phenotype for all 18 traits
traits <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_age_age2_sex_batch_PCs_All.txt", header = T)) 

print(paste0("There are ", nrow(traits), " individuals in the phenotype file and ", ncol(traits)-1, " traits"))

##########################################################################
### STEP 2: make matrix (rows are individuals and columns are traits)

# Turn first column (Individuals into rownames)
rownames(traits) <- traits$IID
traits <- traits[, -1]

	
# Make numeric matrix
traits <- mutate_all(traits, function(x) as.numeric(as.character(x)))

# Tranform in matrix
traits_matrix <- as.matrix(traits)

print(paste0("Final dimensions: ", nrow(traits), " (samples; row) x ", ncol(traits), " (traits; column)"))

#################################################
### STEP 3: Calculate the correlation matrix ###

# Calculate the correlation matrix
traits_matrix_cor <- cor(traits_matrix, use = 'pairwise.complete.obs')
rm(traits_matrix)

# Set missing values to 0 to allow svd()
traits_matrix_cor[which(is.na(traits_matrix_cor))] <- 0

#################################################
### STEP 4: Calculate the eigenvalues ###########

# svd() computes the singular-value decomposition of a rectangular matrix, with $d being a vector containing the singular values of the decomposition
traits_matrix_EV <- svd(traits_matrix_cor)$d 


#################################################
### STEP 5: Calculate Neff ######################

# Neff is defined as the number of eigenvalues required to explain 99.5% of the variation from the traits_matrix data 
sum_EV <- 0
count <- 0

while(sum_EV/sum(traits_matrix_EV) < 0.995) {
	count <- count + 1
	sum_EV <- sum_EV + traits_matrix_EV[count]
}

print(paste0("Neff for ", ncol(traits), " continuous traits: " , count))


#################################################
### Save ########################################
fwrite(data.frame(count), "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Neff/Neff_18continuous_traits.txt", col.names = F, row.names = F, quote = F, sep = "\t")

