# install.packages("BiocManager")
# load related packasges
# for GEO datasets
library(GEOquery) # BiocManager::install("GEOquery")
# for data multipulating
library(tidyverse) # install.packages("tidyverse")
library(dplyr) # install.packages("dplyr")
# for genes' symbols
library(org.Hs.eg.db) # BiocManager::install("org.Hs.eg.db")
# for STRING PPI
library(STRINGdb) # BiocManager::install("STRINGdb")
library(igraph) # install.packages("igraph")

# set work directory
setwd("E:/ResearchData/GEO/GSE30550") # Please replace it to your own working directory

# download GSE30550 datasets
gse <- getGEO("GSE30550", GSEMatrix = TRUE, destdir = ".")

# get expression matrix
expression_data <- exprs(gse[[1]])
rownames(expression_data) <- fData(gse[[1]])$`Gene_ID`

# get phenotype data and seperate subjects and timepoints
pheno_data <- pData(gse[[1]]) %>% 
  separate(title, into = c("subject", "timepoint"), 
           sep = ",", extra = "merge")

# map genes' ENTREZ IDs to genes' symbols
# get genes' ENTREZ IDs
entrez_ids <- fData(gse[[1]])$`Gene_ID` %>% as.character()
# map genes' ENTREZ IDs to genes' symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
symbols_df <- data.frame(ENTREZID = entrez_ids, SYMBOL = gene_symbols) %>% na.omit()

# initialize STRING object(Human)
string_db <- STRINGdb$new(version = "12.0",
                          species = 9606,
                          score_threshold = 700)
# map genes' symbols to STRING IDs
string_mapped <- string_db$map(symbols_df, "SYMBOL", removeUnmappedRows = TRUE)
# get PPI network
ppi_network_initial <- string_db$get_interactions(string_mapped$STRING_id)
# match ENSP IDs to gene symbols in PPI network
ppi_network <- ppi_network_initial %>%
  dplyr::left_join(
    string_mapped %>% dplyr::select(STRING_id, SYMBOL), 
    by = c("from" = "STRING_id")
  ) %>%
  dplyr::rename(from_symbol = SYMBOL) %>%
  dplyr::left_join(
    string_mapped %>% dplyr::select(STRING_id, SYMBOL), 
    by = c("to" = "STRING_id")
  ) %>%
  dplyr::rename(to_symbol = SYMBOL) %>% 
  dplyr::filter(!is.na(from_symbol) & !is.na(to_symbol)) %>%
  dplyr::select(from_symbol, to_symbol, combined_score, from, to)

# construct network and delete single nodes
# construct edges list
edges <- ppi_network %>%
  dplyr::select(from_symbol, to_symbol) %>%
  dplyr::distinct()
# construct graph without direction
g <- graph_from_data_frame(edges, directed = FALSE)
# calculate connected components
components <- components(g)
# tag isolated nodes
isolated_nodes <- names(components$membership)[components$csize[components$membership] <= 1]

ppi_final <- ppi_network %>%
  dplyr::filter(
    !(from_symbol %in% isolated_nodes) & 
      !(to_symbol %in% isolated_nodes)
  )
ppi_final$combined_score <- ppi_final$combined_score/1000

# get usefule expression matrix
unique_symbols <- unique(c(ppi_final$from_symbol, ppi_final$to_symbol))
unique_ids <- string_mapped[string_mapped$SYMBOL %in% unique_symbols,]$ENTREZID
row_names <- c()
for (i in unique_ids){
  row_names <- c(row_names,string_mapped[string_mapped$ENTREZID==i,]$SYMBOL)
}
expression_final <- expression_data[unique_ids,]
rownames(expression_final) <- row_names
# get subjects' list
subjects <- unique(pheno_data$subject)
# get expression matrix for every subject and every timepoint, save
for (i in 1:length(subjects)){
  timepoints <- subset(pheno_data,subject==subjects[i])$timepoint
  geo_accession <- subset(pheno_data,subject==subjects[i])$geo_accession
  expression_subject <- data.frame(expression_final[,geo_accession])
  colnames(expression_subject) <- timepoints
  write.csv(expression_subject,paste0(subjects[i],".csv"),row.names=TRUE)
}
# get PPI network, save
write.csv(ppi_final[,c("from_symbol","to_symbol","combined_score")],
          "D:\\Github\\GitDB\\CPU\\Bachelor\\data\\GEO\\GSE30550\\PPI.csv",row.names=FALSE)

# get genes' function(description)
gene_function <- fData(gse[[1]])[fData(gse[[1]])$`Gene_ID` %in% unique_ids,][,c("Gene_ID","Description")]
rownames(gene_function) <- gene_function$Gene_ID
gene_function <- gene_function[unique_ids,]
gene_function$SYMBOL <- row_names
gene_function <- gene_function[,c("SYMBOL","Gene_ID","Description")]
colnames(gene_function) <- c("SYMBOL","ENTREZID","Description")
# save
write.csv(gene_function,"D:\\Github\\GitDB\\CPU\\Bachelor\\data\\GEO\\GSE30550\\function.csv",row.names=FALSE)

