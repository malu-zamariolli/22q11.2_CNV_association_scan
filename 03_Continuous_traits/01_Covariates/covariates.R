## Covariates extraction

## Covariates 
### age: comes from phenotype files
### age^2
### batches: comes from file /data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt
### sex: comes from the .fam file /data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/plink/_001_ukb_cal_chr1_v2.fam
### PCA: 1- 40 - comes from file /data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt Columns 26-65

############################################################

library(data.table)
library(dplyr)

############# Load data with batch, sex and PCA info #################

# Sample QC
sample_qc <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/geno/ukb_sqc_v2.txt", header = F, select = c(4,26:65), col.names = c("batch", paste0("PC", seq(1,40)))))

# Sample eid
sample_eid <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/plink/_001_ukb_cal_chr1_v2.fam", header = F, select = c(1,5), col.names = c("eid", "sex")))

# Merge
cov <- cbind(sample_eid, sample_qc)


################# Age and age^2 #################

# Identify column containing age
header <- fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb28603.csv", header = T, nrow = 0) # The file with the information is identified with: whereis 21003 in the command line
col_age <- grep("21003-", names(header))

# Extract age and determine age^2
age <- as.data.frame(fread("/data/FAC/FBM/DBC/zkutalik/default_sensitive/UKBB/pheno/ukb28603.csv", header = T, select = c(1, col_age[1]), col.names = c("eid", "age")))
age$age2 <- age$age^2

# Add to main covariate file

cov <- right_join(age, cov, by = "eid") # right_join returns all rows from y (cov) and all columns from x and y (age and cov)

########## Save covariates #############

fwrite(cov, "/scratch/beegfs/FAC/FBM/DBC/zkutalik/default_sensitive/mzamario/Phenotype_extraction/data_input/covariates.txt", col.names = T, row.names = F, quote = F, sep = "\t")



