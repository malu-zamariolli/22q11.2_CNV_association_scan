

####### Steps of the phenotype extraction
# 1. Load file with phenotypes of interest
# 2. List most recent raw UKBB file to each pheno of interest from where phenotypic information will be extracted
# 3. Load non-redacted eids
# 4. Select samples: from file with unrelated individuals 
# 6. Summary info for phenotypes, such as mean, sd and number of NA 
# 7. Number of unique values 
# 8. Inverse Normal transformation for all phenotypes
# 9. Regress out covariates


### Libraries #################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(matrixStats)


#################################################
### Load Main Data #############################

# Open file with the continuous traits to be extracted

pheno_list <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/data_input/continuous_traits_22q112.txt", header = T))

#################################################
### Order raw phenotype files ###################

raw_pheno_file <- file.info(list.files("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno", pattern = "^ukb", full.names = T, recursive = F))

raw_pheno_file <- data.frame(File = rownames(raw_pheno_file), Date = raw_pheno_file[ ,4])

raw_pheno_csv <- raw_pheno_file %>% filter(str_detect(File, "\\.csv"))

raw_pheno_csv <- raw_pheno_csv[order(raw_pheno_csv$Date, decreasing = T), ]
print(paste0("There are ", nrow(raw_pheno_csv), " raw phenotype files"))


# This will list the files from UKBB (with path) and date. It is used later on the script to select the most recent version of a ukbb raw file

#################################################
### Identify most recent location ###############

pheno_list$File <- NA
pheno_list$Col <- NA

for (p in 1:nrow(pheno_list)) {

	# Define pheno
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	print(paste0("Extracting most recent file for ", pheno))
	counter <- 1

	# Define most recent location --> loop through raw files
	while (is.na(pheno_list[p, "File"])) {
	
		header <- fread(as.character(raw_pheno_csv[counter, "File"]), nrow = 0)
		col <- grep(paste0("^", ID, "-"), names(header))
		if (length(col) > 0) {
			pheno_list[p, "File"] <- as.character(raw_pheno_csv[counter, "File"])
			print(paste0("Most recent location: ", sub(".*/", "", pheno_list[p, "File"])))
			pheno_list[p, "Col"] <- paste(col, collapse = "_")}
		counter <- counter +1
	}
} 
rm(p, pheno, ID, counter, header, col)


fwrite(pheno_list, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_list.txt", col.names = T, row.names = F, quote = F, sep = "\t")

##############################################
### Extract phenotypes + Average ################

phenotypes <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/data_input/eid_no_redacted.txt", header = T))

for (p in 1:nrow(pheno_list)) {
	
	# Define pheno
	pheno <- pheno_list[p, "Pheno"] 
	ID <- pheno_list[p, "FieldID"] 
	file <- pheno_list[p, "File"] 
	col <- pheno_list[p, "Col"] 
	print(paste0("Extracting ", pheno, " from ", sub(".*/", "", file)))

	# Extract pheno columns
	temp_pheno <- as.data.frame(fread(file, header = T, select = c(1, as.numeric(str_split(col, pattern = "_")[[1]]))))

	# Average (if more than one column present)
	if (ncol(temp_pheno) > 2) {
	temp_pheno <- data.frame(eid = temp_pheno[, 1], pheno = rowMeans(temp_pheno[, 2:ncol(temp_pheno)], na.rm = T))}

	# Merge
	colnames(temp_pheno) <- c("eid", pheno)
	phenotypes <- full_join(phenotypes, temp_pheno, by = "eid")

}
rm(pheno, ID, file, col, temp_pheno)
print(paste0("Dimensions of phenotype table: ", ncol(phenotypes), " pheno x ", nrow(phenotypes), " eids"))


#################################################
### Replace non-conclusive entries by NA ########

# -1: do not know; -2: only had twins; -3: prefer not to answer
phenotypes[phenotypes < 0] <- NA  # Replaces values smaller than 0 by NA in the entire data.frame, not in a single column

#################################################
### Select samples ##############

# All - Both Sexes
eid_all <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_All.txt", header = F))
print(paste0("There are ", nrow(eid_all), " individuals"))
colnames(eid_all) <- "eid"

pheno_all <- phenotypes[phenotypes$eid %in% eid_all$eid, ]

colnames(pheno_all)[1] <- "IID"
print(paste0("Dimensions of phenotype table for all unrelated individuals: ", ncol(pheno_all), " pheno x ", nrow(pheno_all), " eids"))
fwrite(pheno_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

#################################################
### Summary #####################################

# All
summary_pheno_all <- data.frame(Pheno = colnames(pheno_all)[2:ncol(pheno_all)],
								Num = colSums(!is.na(pheno_all[, c(2:ncol(pheno_all))])),
								NumNA = colSums(is.na(pheno_all[, c(2:ncol(pheno_all))])),
								Mean = colMeans(pheno_all[, c(2:ncol(pheno_all))], na.rm = T),
								Median = colMedians(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								SD = colSds(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Min = colMins(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T),
								Max = colMaxs(as.matrix(pheno_all[, c(2:ncol(pheno_all))]), na.rm = T))
fwrite(summary_pheno_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/summaries/summary_pheno_continuous_raw_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

print("Summary all: done!")
########################################################

print("Number of unique values for each continuous trait")

for (i in 2:ncol(pheno_all)) {
	trait <- colnames(pheno_all)[i]
	print(paste0(trait, " - unique values: ", length(unique(pheno_all[, i]))))
}; rm (i, trait)

print("Starting INT ....")

#################################################
### Inverse Normal Transformation ###############

#  All
pheno_int_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
	pheno <- colnames(pheno_all)[p]
	print(paste0("INT - All: ", pheno))
	pheno_int_all[, pheno] <- qnorm((rank(pheno_all[, p], na.last = "keep")-0.5)/sum(!is.na(pheno_all[, p])))
}; rm (p, pheno)

fwrite(pheno_int_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

#################################################
### Correct for covariates ######################

# Read in covariates
cov <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/data_input/covariates.txt", header = T))
colnames(cov)[1] <- "IID"

# Select subsets
cov_all <- cov[cov$IID %in% eid_all$eid, ] 
print(paste0("Covariates all: ", ncol(cov_all)-1))
print(colnames(cov_all[-1]))


# Regress out covariates - All
pheno_int_cor_all <- data.frame(IID = pheno_int_all$IID)
for (p in 2:ncol(pheno_int_all)) {
	# Define phenotype
	pheno <- colnames(pheno_int_all)[p]
	print(paste0("Covariates - All: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_all, pheno_int_all[, c(1,p)], by = "IID") # returns all rows from y and all columns from x and y
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_int_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_int_cor_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

#################################################
### Correct for covariates on non-INT data ######

# Regress out covariates - All
pheno_cor_all <- data.frame(IID = pheno_all$IID)
for (p in 2:ncol(pheno_all)) {
	# Define phenotype
	pheno <- colnames(pheno_all)[p]
	print(paste0("Covariates (no INT) - All: ", pheno))
	# Create a temporary dataframe
	temp <- right_join(cov_all, pheno_all[, c(1,p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "pheno"
	# Regress out covariates
	pheno_cor_all[, pheno] <- residuals(lm(pheno ~ . , data = temp[, -c(1)], na.action = na.exclude))
}; rm (p, pheno, temp)

fwrite(pheno_cor_all, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_nonINT_age_age2_sex_batch_PCs_All.txt", col.names = T, row.names = F, quote = F, sep = "\t", na = NA)







