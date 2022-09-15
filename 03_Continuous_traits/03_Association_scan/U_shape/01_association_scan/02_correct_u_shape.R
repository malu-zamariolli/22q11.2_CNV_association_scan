# Harmonize GWAS effects to A1 = "T"

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

#################################################
### STEP 1: List files to correct ###############

# These are the output files from the GWAS
list_files <- list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/u_shape/intermediate_files/", pattern = ".glm.linear", full.names = T, recursive = F)
print(paste0("Number of files to be analyzed (number of continuous phenotypes tested): ", length(list_files)))


#################################################
### STEP 2: Loop over files & correct ###########

# Loop over GWAS files
for (f in list_files) {

	# Load the file
	df <- as.data.frame(fread(f, header = T)) 

	# Correct the BETA, L95, U95, T_stat. Hetonly option gives effect of heterozygous (AT) which is Copy-neutral. Change to obtain effect of CNVs
	df[, c(9, 11, 12, 13)] <- -1 * df[, c(9, 11, 12, 13)]
	colnames(df)[c(11,12)] <- c("U95", "L95")
	df <- df[, c(1:10, 12, 11, 13:15)] 
	
	# Save corrected version under the same name for follow-up analyses
	fwrite(df, f, col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

}
