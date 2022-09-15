# Perform a step-wise conditional analysis for the continuous phenotypes
# The SCA will determine the number of independent signals per trait

#################################################
### Libraries ###################################
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(gtools)

#####################################################
####### Create input file for SCA with ##############
### significant hits in files separated by pheno #####

# Define threshold for significance (FWER 5%). 
thr <- 0.05/(6*(113+16))

#Load probes of interest
probes_to_keep <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Probes/504_probes_r2_999.txt", header = F, col.names = "rs"))
print(paste0("There are ", nrow(probes_to_keep), " probes to keep - probes with at least 10 LD-friend (r^2 0.999) "))

# List phenotypes ####

# These are the corrected output files from the GWAS
list_files <- list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/", pattern = ".glm.linear", full.names = T, recursive = F)
pheno_list <- unique(gsub(".*\\.", "", sub(".glm.linear", "", list_files)))
print(paste0("There are ", length(pheno_list), " unique phenotypes in "))
rm(list_files)

# Create an empty dataframe
sig_all <- data.frame()

# Loop over phenotypes
for (p in pheno_list) {

	# Define the phenotype 
	print(paste0("Starting analyzing ", p))

	# List files, merge, and save
	merged_gwas <- data.frame()
	files <- mixedsort(list.files("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/", pattern = paste0("\\.", p, ".glm.linear"), full.names = T, recursive = F))
	print(paste0("Number of files to merge: ", length(files)))
	for (f in files) {
		df <- as.data.frame(fread(f, header = T, select = c(1:6, 8:15)), stringsAsFactors = F)
		colnames(df)[1] <- "CHR"
		merged_gwas <- rbind(merged_gwas, df)}
	
	# Select and save significant hits
	sig <- merged_gwas[merged_gwas$ID %in% probes_to_keep$rs, ]
	sig <- sig[which(sig$P <= thr), ]
	if (nrow(sig) >= 1) {
		sig$PHENO <- p
		sig_all <- rbind(sig_all, sig)
		fwrite(sig, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/gwas_results/SCA_results/significant_", p, "_mirror_after_region_definition.txt.gz"), col.names = T, row.names = T, quote = F, sep = "\t")
		print(paste0(nrow(sig), " significant associations for ", p))
	} else {
		print("No significant association")}

}

print(paste0("Total number of significant associations in : ", nrow(sig_all)))
fwrite(sig_all, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/gwas_results/SCA_results/00_significant_mirror_after_region_definition.txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t")
rm(p, f, files, merged_gwas, sig, sig_all)

#################################################
### SET UP WHILE LOOP ###########################

# Set a parameter to count the number of iterations
round <- 1

# Load and count the number of significant hits from the GWAS
n <- nrow(as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/gwas_results/SCA_results/00_significant_mirror_after_region_definition.txt.gz", header = T, select = c(1:3, 8:9, 13:15))))

# Start a while loop
while (n > 0) {
	print(paste0("STARTING ROUND ", round))

#################################################
### STEP 1: Select top probes ###################

# Access the significant data; For round = 1 the GWAS output; Else the output of the previous SCA iteration
if (round == 1) {
		new_sig <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/gwas_results/SCA_results/00_significant_mirror_after_region_definition.txt.gz", header = T, select = c(1:3, 8:9, 13:15)))
	} else {
		new_sig <- as.data.frame(fread(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/significant_associations/sig_assoc_mirror_SCA_round_after_region_definition", round-1, ".txt.gz"), header = T, select = c(1:3, 8:9, 13:15)))}

# Extract top GW-significant variant/phenotype
new_sig <- new_sig[order(new_sig$P), ]
new_sig <- new_sig[!duplicated(new_sig[, "PHENO"]), ]
new_sig$ROUND <- round
print(paste0("Number of phenotypes with at least one significant association added at round ", round,": ", nrow(new_sig)))


#################################################
### STEP 2: Fill summary file ###################

# Create (round = 1) or load (round > 1) and file the summary file
if (round == 1) {
		sum_file <- new_sig
	} else {
		sum_file <- as.data.frame(fread(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/summary/summary_mirror_SCA_round_after_region_definition", round-1, ".txt"), header = T))
		sum_file <- rbind(sum_file, new_sig)}

print(paste0("Total number of independent signals after round ",round,": ", nrow(sum_file)))
fwrite(sum_file, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/summary/summary_mirror_SCA_round_after_region_definition", round, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)


#################################################
### STEP 3: Probe CNV profile ###################

# Detect the probes whose CNV profile is to be extracted
probes <- sum_file[order(sum_file$ROUND), c(1,3,8,9)]
probes <- probes[!duplicated(probes[, "ID"]), ]
probes <- probes[which(probes$ROUND == round), ]
print(paste0("Number of probes to prepare: ", nrow(probes)))

# Loop over probes to generate their CNV profile
if(nrow(probes) > 0){
for(i in 1:nrow(probes)) {

	# Define the probe
	rs <- as.character(probes[i, "ID"])
	chr <- as.character(probes[i, "CHR"])
	print(paste0("Starting probe: ", rs, " (chr", chr, ")"))

	# Write file
	fwrite(data.frame(rs), paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/cnv_profiles/", rs), col.names = F, row.names = F, quote = F, sep = "\t")

	# Extract profile with PLINK
	input <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/plink_input/mirror/ukb_cnv_chr22") # PLINK file set for the mirror-only momirror
	output <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/cnv_profiles/", rs)
	system(paste("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/softwares/plink/plink",
				 "--bfile", input, 
				 "--keep /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_plink1.9_All.txt",
				 "--extract", output,
				 "--recode tab --out", output))

	# Correct the .ped in R
 	ped <- as.data.frame(fread(paste0(output, ".ped"), header = F, select = c(2,7), col.names = c("IID", "CNV")))
	ped[which(ped$CNV == "A A"), "CNV"] <- -1
	ped[which(ped$CNV == "A T" | ped$CNV == "T A"), "CNV"] <- 0
	ped[which(ped$CNV == "T T"), "CNV"] <- 1
	colnames(ped)[2] <- rs
	fwrite(ped, paste0(output, "_profile"), col.names = T, row.names = F, quote = F, sep = "\t")

	# mirrorete the rs file and the .ped/.map
	unlink(output); unlink(paste0(output, ".log")); unlink(paste0(output, ".map")); unlink(paste0(output, ".ped"))

}
}
rm(i, rs, chr, input, output, ped)


#################################################
### STEP 4: Prepare phenotypes ##################
print("Starting to regress out the covariates")

# Read in INT phenotypes (before covariates correction) for all and select those with significant association; from the phenotype extraction "continuous_trait_extraction.R"
pheno <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/Traits_22q11.2/continuous/pheno_continuous_INT_All.txt", header = T))
pheno <- pheno[ , colnames(pheno) %in% c("IID", unique(new_sig$PHENO))]

# Read in covariates; from the phenotype extraction "continuous_trait_extraction.R"
cov <- as.data.frame(fread("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/data_input/covariates.txt", header = T))
colnames(cov)[1] <- "IID"
cov <- cov[cov$IID %in% pheno$IID, ]

# Correct the phenotypes for covariates
pheno_cor <- data.frame(IID = pheno$IID)
for (i in 2:ncol(pheno)) {

	# Define phenotype and variants to correct for
	p <- colnames(pheno)[i]
	rs_list <- sum_file[which(sum_file$PHENO == p), ]

	# Add SCA variants to the covariate file
	temp <- cov
	for (rs in rs_list$ID) {
		rs_profile <- as.data.frame(fread(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/cnv_profiles/", rs, "_profile"), header = T))
		temp <- left_join(temp, rs_profile, by = "IID")}

	# Add the phenotype of interest to the temp file
	temp <- right_join(temp, pheno[ , names(pheno) %in% c("IID", p)], by = "IID")
	colnames(temp)[ncol(temp)] <- "PHENO"

	# Regress out covariates
	pheno_cor[, p] <- residuals(lm(PHENO ~ . , data = temp[, -c(1)], na.action = na.exclude))
	
}

fwrite(pheno_cor, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/phenotypes/pheno_cor_round_after_region_definition", round, ".txt"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
rm(i, p, rs_list, rs, rs_profile, temp, pheno_cor, new_sig, sum_file)


#################################################
### STEP 5: GWAS ################################


print("Starting GWAS")
system(paste0("mkdir /scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/associations/All/round_after_region_definition", round))

# Define variables
input <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/plink_input/mirror/ukb_cnv_chr22") # PLINK binary file set for mirror-only momirror
pheno_file <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/phenotypes/pheno_cor_round_after_region_definition", round, ".txt") # INT, covariate-corrected phenotype file generate above
sample_file <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/CNV_GWAS/samples/samples_All.txt") # Unrelated  individuals from sample_filtering.R

probe_file <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/After_22q112_definition/Probes/504_probes_r2_999.txt")

output <- paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/associations/All/round_after_region_definition", round,"/gwas_mirror_SCA_round_after_region_definition", round,"_chr22")

chr <- 22
start <- 18630000 
end <- 21910000
  

# GWAS in PLINK v2
	system(paste("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/cauwerx/softwares/plink2/plink2",  	
				  "--bfile", input, 
				  "--pheno",  pheno_file, "--no-psam-pheno", 
				  "--glm omit-ref no-x-sex hide-covar allow-no-covars --ci 0.95",
				  "--keep", sample_file,
				  "--extract", probe_file, 
				  "--out", output))
	

rm(input, pheno_file, sample_file, probe_file, output)


#################################################
### STEP 6: Correct GWAS & ID significants ######

# Threshold for significance (FWER 5%)
thr <- 0.05/(6*(113+16))

# List files to corrrect (harmonization of effect allele to "T")
list_gwas <- list.files(paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/associations/All/round_after_region_definition", round), pattern = ".glm.linear", full.names = T, recursive = F)

# Create an empty dataframe
sig_gwas <- data.frame()

# Loop over files to correct the effect allele and extract significant signals
for (f in list_gwas) {

	# Load the file
	df <- as.data.frame(fread(f, header = T))
	p <- sub(".*\\.", "",sub(".glm.linear", "", f)) 

	# Correct the BETA, CI, A1, REF, and ALT
	df[which(df$A1 == "A"), c(9, 11, 12, 13)] <- -1 * df[which(df$A1 == "A"), c(9, 11, 12, 13)]
	df[which(df$A1 == "A"), "REF"] <- "A"
	df[which(df$A1 == "A"), "ALT"] <- "T"
	df[which(df$A1 == "A"), "A1"] <- "T"
	colnames(df)[c(11,12)] <- c("U95", "L95")
	df <- df[, c(1:10, 12, 11, 13:15)]
	
	# Save corrected version under the same name for follow-up analyses
	fwrite(df, f, col.names = T, row.names = F, quote = F, sep = "\t", na = NA)

	# Extract significant signals
	df <- df[, c(1:6, 8:15)]
	colnames(df)[1] <- "CHR"
	df$PHENO <- p
	temp_sig <- df[which(df$P <= thr), ]
	sig_gwas <- rbind(sig_gwas, temp_sig)

}

print(paste0("Number of significant associations after round ", round, ": ", nrow(sig_gwas)))
fwrite(sig_gwas, paste0("/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Traits_22q112_region_CNV_GWAS/Continuous/mirror/intermediate_files/significant_associations/sig_assoc_mirror_SCA_round_after_region_definition", round, ".txt.gz"), col.names = T, row.names = F, quote = F, sep = "\t", na = NA)
rm(f, df, temp_sig)


#################################################
### STEP 7: Update loop variables ###############

round <- round + 1
n <- nrow(sig_gwas)
rm(sig_gwas)

}
