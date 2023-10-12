## PREPROCESSING SCRIPT FOR CANINE AND HUMAN PCA PLOT ## 
## PRIMARY PURPOSE TO GET MASTER LOGGED EXPRESSION MATRICES FOR CANINE AND HUMAN
## AND A MASTER COMBINED METADATA FILE WITH SIMPLIFIED CANCER NAMES ##  
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(here)

#### CANINE PREPROCESSING #### 
#get the project ids and master metadata file 
projects <- read_table(here("data","canis_familiaris","BioProjIDs.txt"))
excluded <- c("PRJNA914497","PRJNA338759") 
projects <- projects[!(projects$BioProjIDs %in% excluded),]
canine_metadata <- readr::read_tsv(here("data","canis_familiaris","metadata_final.txt"))

# create expression matrix # 
count_files <- here("data","canis_familiaris",projects$BioProjIDs,"freya_results","dexseq_count",paste0(projects$BioProjIDs,".count"))
all(file.exists(count_files))
canine_expr <- lapply(count_files, readr::read_tsv)
canine_expression <- canine_expr %>% purrr::reduce(inner_join)
dim(canine_expression)

# correct hugo names # 
library(biomaRt)
old_ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", host = "https://may2021.archive.ensembl.org",
                          dataset = "clfamiliaris_gene_ensembl",version = 104)
gene_names <-  getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = list(ensembl_gene_id = canine_expression$Genes),
                     mart = old_ensembl)
canine_final <- merge(gene_names, canine_expression, by.x = "ensembl_gene_id", by.y = "Genes")
canine_final <- canine_final[,-1] #drop entrez ids 
canine_final[1:5,1:5]
write.table(canine_final, here("data","testing","canine_expr.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)

# summarise values by gene (for duplicate genes) # 
canine_final <- canine_final %>% group_by(external_gene_name) %>% summarise_all(mean)
canine_final <- column_to_rownames(canine_final, "external_gene_name")
canine_final <- canine_final[-1,] # first row was gene not matched to symbol, so got rid of it

# log gene values # 
canine_final <- log2(canine_final + 1)
canine_final[1:5,1:5]
canine_final <- rownames_to_column(canine_final, "Gene")
canine_final[1:5,1:5]
write.table(canine_final,here("data","processed_data","canine_logged_expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

# create batch ids # 
batch.ids <- unlist(sapply(seq_along(projects$BioProjIDs), function(i) {rep(projects$BioProjIDs[i],ncol(canine_expr[[i]])-1)} ))
temp <- colnames(canine_expression)
head(temp)
temp <- temp[-1]
names(batch.ids) <- temp

#write batch.ids to file 
batch.dat <- data.frame(batch.ids,names(batch.ids))
head(batch.dat)
colnames(batch.dat) <- c("Batch", "Sample")
write.table(batch.dat, here("data","processed_data","canine_batch_ids.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

#add batch ids to metadata file 
temp_canine_met <- merge(batch.dat, canine_metadata, by.x = "Sample", by.y = "Run")
write.table(temp_canine_met,here("data","processed_data","final_canine_metadata.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

#### HUMAN PAN-CANCER PREPROCESSING ####
#get list of all pancancer diseases 
all_diseases_init <- list.files(here("data","homo_sapiens","panCancer_Atlas"))
all_diseases <- all_diseases_init[grep("2018",all_diseases_init)]

# get metadata 
human_metadata <- readr::read_tsv(here("data","homo_sapiens","panCancer_Atlas",all_diseases,"data_clinical_sample.txt"))
# the column names are weird, we actually want the 4th row read in to be the column names 
correct_cols <- human_metadata[4,]
colnames(human_metadata) <- correct_cols

# get pre made pancancer expression matrix - made with script: merge_pancancer_data.R
human_expression <- readr::read_tsv(here("data","homo_sapiens","panCancer_Atlas","final_all_panCancerExpression.tsv"))
dim(human_metadata)
human_metadata[1:5,1:5]
human_expression[1:5,1:5]

# get rid of the exrtra note lines that were at the top of the metadata in each cancer file 
human_metadata <- human_metadata[human_metadata$SAMPLE_ID %in% colnames(human_expression), ]
dim(human_metadata)

#create batch ids 
human_expr <- here("data","homo_sapiens","panCancer_Atlas",all_diseases,"data_mrna_seq_v2_rsem.txt")
all(file.exists(human_expr))
human_expre <- lapply(human_expr,readr::read_tsv)

#create human batch ids - this isnt pretty but I couldnt figure out a better way
# (p.s. we substract 2 this time becuase we have both hugo and entrz names columns)
human_batch.ids <- unlist(sapply(seq_along(all_diseases), function(i) {rep(all_diseases[i],ncol(human_expre[[i]]) - 2)} ))
length(human_batch.ids)
sample_id <- colnames(human_expression)
sample_id <- sample_id[-1] #gets rid of entrez names 
head(sample_id)
names(human_batch.ids) <-sample_id
human_batch_dat <- data.frame(human_batch.ids, names(human_batch.ids))
head(human_batch_dat)
colnames(human_batch_dat) <- c("Batch","Sample")
write.table(human_batch_dat, here("data", "processed_data","human_batch_ids.tsv"), sep="\t", col.names = TRUE, row.names = FALSE)

# add batch ids to metadata in case helpful
temp_human_met <- merge(human_batch_dat, human_metadata, by.x = "Sample", by.y = "SAMPLE_ID")
write.table(temp_human_met,here("data","processed_data","final_human_metadata.txt"), sep = "\t", col.names = TRUE, row.names = FALSE)

# log expression and save it
human_expression <- column_to_rownames(human_expression, "Hugo_Symbol")
human_expression <- log2(human_expression + 1)
human_expression[1:5,1:5]
#replace NA with 0 
human_expression[is.na(human_expression)] <- 0 
human_expression <- rownames_to_column(human_expression, "Gene")
human_expression[1:5,1:5]
write.table(human_expression,here("data","processed_data","human_logged_expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

## COMBINE METADATA ## 
head(temp_canine_met)
canine_cut <- temp_canine_met %>% dplyr::select(Sample, Batch, BioProject, Tumor, CancerType)
human_cut <- temp_human_met %>% dplyr::select(Sample, Batch, ONCOTREE_CODE, TUMOR_TYPE,CANCER_TYPE)
colnames(canine_cut) <- c("Sample", "Batch", "ProjCode","Tumor","CancerType")
colnames(human_cut) <- c("Sample", "Batch", "ProjCode","Tumor","CancerType")
canine_cut$Species <- "Canine"
human_cut$Species <- "Human"

all_meta <- rbind(human_cut,canine_cut)
dim(all_meta)

# correct cancer names 
table(all_meta$CancerType)

# - seperated by groups 
## lymphoma - theses are all canine samples 
all_meta$CancerType[grep("T-cell Lymphoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Lymphoma"
all_meta$CancerType[grep("B-cell Lymphoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Lymphoma"
all_meta$CancerType[grep("B-cell Lumphoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Lymphoma"

## bladder cancer, changed canine label to match human
all_meta$CancerType[grep("Bladder Tumor", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Bladder Cancer"

## breast cancer changed canine to match human 
all_meta$CancerType[grep("Mammary Tumor", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Breast Cancer"

## brain cancer - created brain cancer column, meningioma unique ot canine
all_meta$CancerType[grep("Glioblastoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Brain Cancer"
all_meta$CancerType[grep("Glioma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Brain Cancer"
all_meta$CancerType[grep("Meningioma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Brain Cancer"

## Melanoma
all_meta$CancerType[grep("Cutaneous Melanoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Melanoma"
all_meta$CancerType[grep("Acral Melanoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Melanoma"
all_meta$CancerType[grep("Oral Melanoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Melanoma"
## ASK KILEY ABOUT BELOW GROUPING 
all_meta$CancerType[grep("Mucosal Melanoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Melanoma"
all_meta$CancerType[grep("Ocular Melanoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Melanoma"

#sarcoma - soft tissue sarcoma is only cainice and Hemangiosarcoma is only canince and osteosarcoma only canine 
all_meta$CancerType[grep("Soft Tissue Sarcoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Sarcoma"
all_meta$CancerType[grep("Osteosarcoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Sarcoma"
all_meta$CancerType[grep("Hemangiosarcoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Sarcoma"

#carcinoma 
all_meta$CancerType[grep("Renal Non-Clear Cell Carcinoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Renal Carcinoma"
all_meta$CancerType[grep("Renal Clear Cell Carcinoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Renal Carcinoma"
#all_meta$CancerType[grep("Cholangiocarcinoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Carcinoma"

#head and neck cancer - changed to match human -CHECK WITH KILEY COULD ALSO MAKE ALL THEESE Head and Neck Squamous Cell Carcinoma  
all_meta$CancerType[grep("Head and Neck Squamous Carcinoma", all_meta$CancerType, ignore.case = FALSE, fixed = TRUE)] <- "Head and Neck Cancer"

write.table(all_meta, here("data","processed_data","canine_human_merged_metadata.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)


