#Author: Kaido Lepik

library(readxl)
library(tidyverse)
library(data.table)
library(RSelenium)
library(Rcpp)

setwd("BEPE_UNIL/Project/Traits_22q11.2/Web_scraping/")

#############################################################
### Open df with continuous traits

dat <- as.data.frame(fread("continuous_fieldID_JURA.txt", na.strings = ""))

# Filter traits that are on JURA
dat <- dat %>% filter(!(is.na(File))) # 1769 rows

#################################################################
## Web scraping continuous traits

scrape_HPO_element_text <- function(driver, selector) {
    HPO_text <- driver$findElements("css selector", selector) %>%
        lapply(function(x) x$getElementText()) %>%
        unlist()
    
    return(HPO_text)
}

scrape_HPO_matches <- function(driver, UKBB_label, label_type = c("term", "disease"), identifier_selector = c("mat-cell.mat-column-ontologyId", "mat-cell.mat-column-diseaseId"), 
                               name_selector = c("mat-cell.mat-column-name", "mat-cell.mat-column-dbName"), sleeping_time = 1) {
    # Make sure that function arguments are correct
    label_type <- match.arg(label_type)
    identifier_selector <- match.arg(identifier_selector)
    name_selector <- match.arg(name_selector)
    
    # Navigate to the URL and give it some time to load the page
    url_string <- paste0("https://hpo.jax.org/app/browse/search?q=", str_replace_all(UKBB_label, " ", "%20"), "&navFilter=", label_type)
    driver$navigate(url_string)
    Sys.sleep(sleeping_time)
    
    # Gather IDs and names for all the HPO label matches
    HPO_identifiers <- scrape_HPO_element_text(driver, identifier_selector)
    HPO_names <- scrape_HPO_element_text(driver, name_selector)
    
    # Return the results as a tibble
    tibble(UKBB_label = UKBB_label, HPO_label_type = label_type, HPO_identifier = HPO_identifiers, HPO_name = HPO_names)
}

HPO_matches_to_UKBB <- function(driver, UKBB_label) {
    term_matches <- scrape_HPO_matches(driver, UKBB_label, label_type = "term", identifier_selector = "mat-cell.mat-column-ontologyId", name_selector = "mat-cell.mat-column-name")
    disease_matches <- scrape_HPO_matches(driver, UKBB_label, label_type = "disease", identifier_selector = "mat-cell.mat-column-diseaseId", name_selector = "mat-cell.mat-column-dbName")
    
    rbindlist(list(term_matches, disease_matches), fill = TRUE)
}


# Start the Selenium server and remote browser
rD <- rsDriver(browser = "chrome", chromever = "97.0.4692.71") 
driver <- rD[["client"]]

# Results in a list
HPO_matches_continuousUKBB <- lapply(dat$Pheno, function(UKBB_label) HPO_matches_to_UKBB(driver, UKBB_label)) %>%
    set_names(dat$Pheno)

# Results in a data table
dat_HPO_matches_continuousUKBB <- rbindlist(HPO_matches_continuousUKBB, fill = TRUE)

write.table(dat_HPO_matches_continuousUKBB, "HPO_continuousUKB_webscraping_all_result.txt", sep = "\t", quote = F,
            col.names = T, row.names = F)

# Filter the ones that obtained results

dat_final <- dat_HPO_matches_continuousUKBB %>% filter(!is.na(HPO_identifier))

write.table(dat_final, "HPO_continuousUKB_webscraping_result.txt", sep = "\t", quote = F,
            col.names = T, row.names = F)

# Save list of those that returned without HPO matches in web scraping
dat_no <- dat_HPO_matches_continuousUKBB %>% filter(is.na(HPO_identifier))
dat_no <- unique(dat_no[, c(1,3,4)])

write.table(dat_no, "HPO_continuousUKB_webscraping_results_no_match.txt", sep = "\t", quote = F,
            col.names = T, row.names = F)

rD$server$process$kill() # Kill the server and the browser (and free up the used port)

#######################################################################################

### Step 2

data <- read_xlsx("HPO_continuousUKB_webscraping_results_no_match_manual.xlsx")
dat <- data.frame(keyword = unique(data$keyword))

## Step 2: Web scraping continuous traits

scrape_HPO_element_text <- function(driver, selector) {
    HPO_text <- driver$findElements("css selector", selector) %>%
        lapply(function(x) x$getElementText()) %>%
        unlist()
    
    return(HPO_text)
}

scrape_HPO_matches <- function(driver, UKBB_label, label_type = c("term", "disease"), identifier_selector = c("mat-cell.mat-column-ontologyId", "mat-cell.mat-column-diseaseId"), 
                               name_selector = c("mat-cell.mat-column-name", "mat-cell.mat-column-dbName"), sleeping_time = 1) {
    # Make sure that function arguments are correct
    label_type <- match.arg(label_type)
    identifier_selector <- match.arg(identifier_selector)
    name_selector <- match.arg(name_selector)
    
    # Navigate to the URL and give it some time to load the page
    url_string <- paste0("https://hpo.jax.org/app/browse/search?q=", str_replace_all(UKBB_label, " ", "%20"), "&navFilter=", label_type)
    driver$navigate(url_string)
    Sys.sleep(sleeping_time)
    
    # Gather IDs and names for all the HPO label matches
    HPO_identifiers <- scrape_HPO_element_text(driver, identifier_selector)
    HPO_names <- scrape_HPO_element_text(driver, name_selector)
    
    # Return the results as a tibble
    tibble(UKBB_label = UKBB_label, HPO_label_type = label_type, HPO_identifier = HPO_identifiers, HPO_name = HPO_names)
}

HPO_matches_to_UKBB <- function(driver, UKBB_label) {
    term_matches <- scrape_HPO_matches(driver, UKBB_label, label_type = "term", identifier_selector = "mat-cell.mat-column-ontologyId", name_selector = "mat-cell.mat-column-name")
    disease_matches <- scrape_HPO_matches(driver, UKBB_label, label_type = "disease", identifier_selector = "mat-cell.mat-column-diseaseId", name_selector = "mat-cell.mat-column-dbName")
    
    rbindlist(list(term_matches, disease_matches), fill = TRUE)
}


# Start the Selenium server and remote browser
rD <- rsDriver(browser = "chrome", chromever = "97.0.4692.71") # Esse aqui está indo


driver <- rD[["client"]]

# Results in a list
HPO_matches_continuousUKBB <- lapply(dat$keyword, function(UKBB_label) HPO_matches_to_UKBB(driver, UKBB_label)) %>%
    set_names(dat$keyword)

# Results in a data table
dat_HPO_matches_continuousUKBB <- rbindlist(HPO_matches_continuousUKBB, fill = TRUE)
colnames(dat_HPO_matches_continuousUKBB)[[1]] <- "keyword"
# Merge with continuous traits
data <- data[, c(1,2)]
dat_result <- merge(dat_HPO_matches_continuousUKBB, data, by = "keyword", allow.cartesian = T)

write.table(dat_result, "HPO_continuousUKB_webscraping_all_result_2.txt", sep = "\t", quote = F,
            col.names = T, row.names = F)

# Filter the ones that obtained results

dat_final <- dat_result %>% filter(!is.na(HPO_identifier))
dat_final <- unique(dat_final)

write.table(dat_final, "HPO_continuousUKB_webscraping_result_2.txt", sep = "\t", quote = F,
            col.names = T, row.names = F)

rD$server$process$kill() # Kill the server and the browser (and free up the used port)

