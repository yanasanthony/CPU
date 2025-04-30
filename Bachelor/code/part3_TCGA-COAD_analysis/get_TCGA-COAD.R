# install.packages("BiocManager")
# load related packasges
# for TCGA datasets
library(TCGAbiolinks) # BiocManager::install("TCGAbiolinks")
# for data multipulating
library(SummarizedExperiment) # BiocManager::install("SummarizedExperiment")
library(dplyr) # install.packages("dplyr")
library(XML) # install.packages("XML")
# for STRING PPI
library(STRINGdb) # BiocManager::install("STRINGdb")
library(igraph) # install.packages("igraph")

# set working directory
setwd("E:\\ResearchData\\GDC") # Please replace it to your own working directory 

###### download TCGA-COAD transcriptome
# query GDC database for TCGA-COAD transcriptome data
query_exp <- GDCquery(
  project="TCGA-COAD",
  data.category="Transcriptome Profiling",
  data.type="Gene Expression Quantification",
  workflow.type="STAR - Counts"
)
# download
GDCdownload(query_exp,method="api",files.per.chunk=50)
# get data
exp_data <- GDCprepare(query_exp)

# extract transcriptome:TPM
TPM <- assay(exp_data,"tpm_unstrand")
keep_genes <- rowSums(TPM>=1) >= 0.5*ncol(TPM)
TPM <- TPM[keep_genes,]
# extract gene info and gene name
gene_info <- rowData(exp_data)
gene_ids <- gene_info$gene_id
gene_symbols <- gene_info$gene_name
symbols_df <- data.frame(ENSEMBLID = gene_ids, SYMBOL = gene_symbols) %>% 
  filter(!grepl("_PAR_Y$", ENSEMBLID)) %>% na.omit()
symbols_df <- symbols_df[symbols_df$ENSEMBLID %in% rownames(TPM),]

# find duplicated symbols
duplicated_df <- symbols_df %>%
  group_by(SYMBOL) %>%
  filter(n() > 1) %>%
  ungroup()
# calculate mean to replace original rows
duplicated_symbols <- unique(duplicated_df$SYMBOL)
for (symbol in duplicated_symbols){
  ids <- duplicated_df[duplicated_df$SYMBOL==symbol,]$ENSEMBLID
  dup <- TPM[ids,]
  TPM <- TPM[!(rownames(TPM) %in% ids),]
  means <- colMeans(dup)
  TPM <- rbind(TPM,means)
  rownames(TPM)[nrow(TPM)] <- ids[1]
}
# update symbols_df
symbols_df <- symbols_df[symbols_df$ENSEMBLID %in% rownames(TPM),]
# update TPM
TPM <- TPM[rownames(TPM) %in% symbols_df$ENSEMBLID,]


## initialize STRING object(Human)
string_db <- STRINGdb$new(version = "12.0",
                          species = 9606,
                          score_threshold = 700)
# map genes' symbols to STRING IDs
string_mapped <- string_db$map(symbols_df, "SYMBOL", removeUnmappedRows = TRUE)
string_mapped <- string_mapped %>%
  distinct(STRING_id,.keep_all=TRUE) %>%
  distinct(ENSEMBLID,.keep_all=TRUE)
# update symbols_df and TPM
symbols_df <- symbols_df[symbols_df$ENSEMBLID %in% string_mapped$ENSEMBLID,]
TPM <- TPM[rownames(TPM) %in% symbols_df$ENSEMBLID,]

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
  ) %>% unique()
ppi_final$combined_score <- ppi_final$combined_score/1000

# get usefule expression matrix
unique_symbols <- unique(c(ppi_final$from_symbol, ppi_final$to_symbol))
unique_ids <- string_mapped[string_mapped$SYMBOL %in% unique_symbols,]$ENSEMBLID
row_names <- c()
for (i in unique_ids){
  row_names <- c(row_names,string_mapped[string_mapped$ENSEMBLID==i,]$SYMBOL)
}
TPM_final <- TPM[unique_ids,]
rownames(TPM_final) <- row_names

# reference(normal) samples and cancer samples, totally 524
samples <- getResults(query_exp,cols=c("cases"))
# for 481 cancer samples
cancer_samples <- TCGAquery_SampleTypes(barcode = samples,
                                        typesample = "TP")
case_TPM <- TPM_final[,cancer_samples]
# for 43 reference(normal) samples
ref_samples <- TCGAquery_SampleTypes(barcode = samples,
                                     typesample = "NT")
ref_TPM <- TPM_final[,ref_samples]



###### download TCGA-COAD clinical information
# query GDC database for TCGA-COAD clinical information
query_clin <- GDCquery(
  project="TCGA-COAD",
  data.category="Clinical",
  data.type="Clinical Supplement",
  data.format="BCR xml")
# download
GDCdownload(query_clin,method="api",files.per.chunk=50)
# get data
# general information(age,sex,race) for 459 patients
patients_data <- GDCprepare_clinic(query_clin,
                                   clinical.info="patient") %>%
  select(where(~!all(is.na(.x)))) %>% arrange(stage_event_pathologic_stage)
# follow up information for 459 patients
follow_up_data <- GDCprepare_clinic(query_clin,
                                    clinical.info="follow_up") %>%
  select(where(~!all(is.na(.x))))
# stage event information for 459 patients
stage_event_data <- GDCprepare_clinic(query_clin,
                                      clinical.info="stage_event") %>%
  select(where(~!all(is.na(.x))))
# drug information for 459 patients
drug_data <- GDCprepare_clinic(query_clin,
                               clinical.info="drug") %>%
  select(where(~!all(is.na(.x))))
# radiation information for 459 patients
radiation_data <- GDCprepare_clinic(query_clin,
                                    clinical.info="radiation") %>%
  select(where(~!all(is.na(.x))))

# conbine baseline information
# get barcode
barcode <- patients_data$bcr_patient_barcode
# get stage
stage <- left_join(
  data.frame(bcr_patient_barcode = barcode),
  stage_event_data,
  by = "bcr_patient_barcode"
) %>% pull(pathologic_stage) %>% as.character() %>% na_if("")
names(stage) <- barcode
# get age
age <- floor(abs(patients_data$days_to_birth)/365.25)
names(age) <- barcode
# get gender
gender <- patients_data$gender %>% as.character()
names(gender) <- barcode
# get race or enthnic group
race <- patients_data$race_list %>% as.character() %>% na_if("")
names(race) <- barcode
enthnicity <- patients_data$ethnicity %>% as.character() %>% na_if("")
names(enthnicity) <- barcode
# get BMI
height <- patients_data$height #cm
weight <- patients_data$weight #g
BMI <- (weight/2)/(height/100)**2
names(BMI) <- barcode
# get histological type
histological_type <- patients_data$histological_type %>% as.character()
names(histological_type) <- barcode
# get drug use
drug <- c()
for (code in barcode) {
  valid_therapy <- c()
  if (!code %in% unique(drug_data$bcr_patient_barcode)){
    valid_therapy <- NA
  }else{
    therapy <- drug_data[drug_data$bcr_patient_barcode==code,]$therapy_types %>% 
      as.character() %>%
      na_if("")
    for (element in therapy) {
      if (!(element %in% c(NA,"Ancillary", "Chemotherapy", "Immunotherapy", "Targeted Molecular therapy"))) {
        valid_therapy <- c(valid_therapy,"Ancillary")
      }else{valid_therapy <- c(valid_therapy,element)}
    }
    valid_therapy <- unique(valid_therapy)
    valid_therapy <- paste(valid_therapy,collapse=",")
  }
  drug <- c(drug,valid_therapy)
}
names(drug) <- barcode
# get radiation use
radiation <- c()
for (code in barcode) {
  if (!code %in% unique(drug_data$bcr_patient_barcode)){
    radiation <- c(radiation,"No")
  }else{
    radiation <- c(radiation,"Yes")
  }
}
names(radiation) <- barcode
# get cancer status
cancer_status <- patients_data$person_neoplasm_cancer_status %>% as.character() %>% na_if("")
names(cancer_status) <- barcode
# get times when event occurred and events
times <- c()
events <- c()
death <- patients_data$days_to_death
follow <- patients_data$days_to_last_followup
for (i in 1:length(death)) {
  if (is.na(death[i]) & is.na(follow[i])) {
    times <- c(times, NA)
    events <- c(events, NA)
  } else if (is.na(death[i]) & !is.na(follow[i])) {
    times <- c(times, follow[i])
    events <- c(events, 0)
  } else if (!is.na(death[i]) & is.na(follow[i])) {
    times <- c(times, death[i])
    events <- c(events, 1)
  } else {
    times <- c(times, min(death[i], follow[i]))
    events <- c(events, ifelse(death[i] <= follow[i], 1, 0))
  }
}

# match samples to timepoints
ranksamples <- c()
for (code in barcode){
  for (i in colnames(case_TPM)){
    uuid <- substr(i,1,12)
    if (uuid == code){
      ranksamples <- c(ranksamples,i)
    }
  }
}
case_TPM <- case_TPM[,ranksamples]

###### save
# save case TPM
write.csv(case_TPM,"D:\\Github\\GitDB\\CPU\\Bachelor\\data\\TCGA\\TCGA-COAD\\case\\case.csv",
          row.names=TRUE)
# save reference TPM
write.csv(ref_TPM,"D:\\Github\\GitDB\\CPU\\Bachelor\\data\\TCGA\\TCGA-COAD\\ref\\ref.csv",
          row.names=TRUE)
# save clinical information
clinical_info <- data.frame(
  stage=stage,
  age=age,
  gender=gender,
  race=race,
  enthnicity=enthnicity,
  BMI=BMI,
  histological_type=histological_type,
  therapy=drug,
  ratiation=radiation,
  event=events,
  time=times)
rownames(clinical_info) <- barcode
write.csv(clinical_info,"D:\\Github\\GitDB\\CPU\\Bachelor\\data\\TCGA\\TCGA-COAD\\clinical\\clinical_info.csv")

# save PPI network
write.csv(ppi_final[,c("from_symbol","to_symbol","combined_score")],
          "D:\\Github\\GitDB\\CPU\\Bachelor\\data\\TCGA\\TCGA-COAD\\PPI.csv",row.names=FALSE)
