# Install related packages, mainly BiocManager and TCGAbiolinks
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("TCGAbiolinks", quietly = TRUE)){
  BiocManager::install("TCGAbiolinks")
}

# load related packasges
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(XML)

# Define a function to download transcriptome data
download_rna_data <- function(download_path, cancer_type, output_dir){
  
  # parameters:
  # download_path: path to GDC cart files downloaded
  # cancer_type: input cancer type you want to do research, eg. COAD (according to TCGA)
  # output_dir: path to csv files
  
  message(paste0("Starting download transcriptome data of TCGA-", cancer_type))
  
  # Set work path
  setwd(download_path)
  
  # Download transcriptome
  query_exp <- GDCquery(
    project = paste0("TCGA-", cancer_type),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query_exp)
  exp_data <- GDCprepare(query_exp)
  
  # Extract transcriptome: TPM, FPKM, FPKM-UQ
  assays_list <- list(
    TPM = assay(exp_data, "tpm_unstrand"),
    FPKM = assay(exp_data, "fpkm_unstrand"),
    FPKM_UQ = assay(exp_data, "fpkm_uq_unstrand")
  )
  
  # Extract gene info and gene name
  gene_info <- rowData(exp_data)
  gene_id <- gene_info$gene_id
  gene_name <- gene_info$gene_name
  
  # Create output directory, if not exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (data_type in names(assays_list)) {
    data <- assays_list[[data_type]]
    colnames(data) <- colnames(data)
    rownames(data) <- gene_id
    data <- cbind(gene_name, data)
    data <- cbind(gene_id = rownames(data), data)
    write.csv(data, 
              file.path(output_dir, paste0(data_type, "_transcriptome.csv")), 
              row.names = FALSE)
  }
  
  message("Transcriptome has been saved in CSV!")
}

# Define a function to download clinical data
download_clin_data <- function(download_path, cancer_type, output_dir){
  
  # parameters:
  # download_path: path to GDC cart files downloaded
  # cancer_type: input cancer type you want to do research, eg. COAD (according to TCGA)
  # output_dir: path to csv files
  
  message(paste0("Starting download clinical data of TCGA-", cancer_type))
  
  # Set work path
  setwd(download_path)
  
  # Download clinical data
  query_clin <- GDCquery(project = paste0("TCGA-", cancer_type),
                         data.category = "Clinical",
                         data.type = "Clinical Supplement",
                         data.format = "BCR xml")
  GDCdownload(query_clin)
  clinical <- GDCquery_clinic(project = paste0("TCGA-", cancer_type), 
                              type = "clinical")
  clinical <- clinical[, -1]
  
  
  # Save clinical in csv
  write.csv(clinical,
            file.path(output_dir, paste0(cancer_type,"_clinical.csv")),
            row.names = FALSE)
  
  message("Clinical data has been saved in CSV!")
}

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure that the correct number of arguments are passed
if(length(args) != 3) {
  stop("Error: Please provide download_path, cancer_type and output_dir as arguments")
}

# Assign arguments to variables
download_path <- args[1]
cancer_type <- args[2]
output_dir <- args[3]

# Call functions
download_rna_data(download_path, cancer_type, output_dir)
download_clin_data(download_path, cancer_type, output_dir)
