# install package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TCGAbiolinks")
browseVignettes("TCGAbiolinks")

# load packages
library(TCGAbiolinks)
library(dplyr)
library(reshape2)
library(stringr)

# settings
work_dir <- getwd()
data_dir <- paste0(work_dir, '/GDC/', gsub("-","_","TCGA-BRCA"))
##################################### microRNA data ##################################### 
# get TCGA barcode and data statistics
query <- GDCquery(project = "TCGA-BRCA",
                  data.category =  "Transcriptome Profiling",
                  data.type = "Isoform Expression Quantification",
                  legacy = FALSE,
                  workflow.type = "BCGSC miRNA Profiling")
samplesDown <- getResults(query, cols = c("cases"))
cat('总样本数量：', length(samplesDown))

table(str_sub(samplesDown,14,15))

dataSmTP <- TCGAquery_SampleTypes(samplesDown, "TP")
cat('总TP样本数量：', length(dataSmTP))

# donwload TCGA microRNA data files
file_data_dir <- paste0(data_dir, "_", "miRNA.rda")
GDCdownload(query, method = "client", directory = data_dir, files.per.chunk = 6)

# combine data files to one file and extract data
raw_data <- GDCprepare(query, save = TRUE, directory = data_dir, save.filename = file_data_dir)

miRNA_id_mature <- read.table('miRNA_mature.txt', header = TRUE, stringsAsFactors = FALSE)
new_data <- raw_data %>% 
#  mutate(miRBase_id = str_match(miRNA_region,'\\w+,(\\w+)')[, 2]) %>%
  mutate(miRBase_id = str_extract(miRNA_region,'[^,]+$')) %>% 
  filter(str_detect(miRBase_id,'^M')) %>% 
  left_join(miRNA_id_mature)

read_countdata <- dcast(new_data, miRNA_mature ~ barcode, sum, value.var = 'read_count') %>% 
  filter(!is.na(miRNA_mature))

# save to local file
write.table(read_countdata, paste0(data_dir, "_", "miRNA_counts", ".txt"), sep = '\t', 
            row.names = F, quote = F, col.names = T)



##################################### Clinical data #####################################
# TCGA clinical data statistics
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  file.type = "xml",
                  legacy = FALSE)
samplesDown <- getResults(query, cols = c("cases"))
cat('总样本数量：', length(samplesDown))

# donwload TCGA clinical data files
file_data_dir <- paste0(data_dir, "_", "clinical.txt")
GDCdownload(query, method = "client", directory = data_dir, files.per.chunk = 6)

# combine data files to one file and extract data
clinical.patients <-  GDCprepare_clinic(query, clinical.info = "patient", directory = data_dir)

# save to local file
write.table(clinical.patients, paste0(data_dir, "_", "clinical.patients", ".txt"), sep = '%', 
            row.names = F, quote = F, col.names = T, na = 'NA')
