# Install related packages, mainly BiocManager and GEOquery
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("GEOquery", quietly = TRUE)){
  BiocManager::install("GEOquery")
}

# load related packasges
library(GEOquery)
library(tidyverse)
library(limma)

# set work directory
setwd("E:/ResearchData/BachelorThesis/GEO/GSE30550")

# download GSE30550 datasets
gse <- getGEO("GSE30550", GSEMatrix = TRUE, destdir = ".")

# get expression matrix
expression_data <- exprs(gse[[1]])
rownames(expression_data) <- fData(gse[[1]])$`Gene_ID`

# get phenotype data
pheno_data <- pData(gse[[1]])
pheno_data <- pData(gse[[1]]) %>% 
  separate(title, into = c("subject", "timepoint"), 
           sep = ",\\s+", extra = "merge")

# get subjects' list
subjects <- unique(pheno_data$subject)
# get expression matrix for every subject and every timepoint, save
for (i in 1:length(subjects)){
  timepoints <- subset(pheno_data,subject==subjects[i])$timepoint
  geo_accession <- subset(pheno_data,subject==subjects[i])$geo_accession
  expression_subject <- data.frame(expression_data[,geo_accession])
  colnames(expression_subject) <- timepoints
  write.csv(expression_subject,paste0(subjects[i],".csv"),row.names = TRUE)
}

