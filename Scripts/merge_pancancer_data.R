# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(here)

all_diseases_init <- list.files(here("data","homo_sapiens","panCancer_Atlas")) 
all_diseases <- all_diseases_init[grep("2018",all_diseases_init)]
human_expr <- here("data","homo_sapiens","panCancer_Atlas",all_diseases,"data_mrna_seq_v2_rsem.txt")
all(file.exists(human_expr))
human_expre <- lapply(human_expr,readr::read_tsv)

# #first get list of genes that is common to all datasets
inital_list <- human_expre[1][[1]]$Hugo_Symbol
inital_list <- inital_list[!(is.na(inital_list))]
for(i in 2:length(human_expre)){
  curr <- human_expre[i][[1]]
  inital_list <- intersect(inital_list, curr$Hugo_Symbol)
}
# 
# #start dat 
all_dat <- human_expre[1][[1]]
all_dat <- all_dat %>% dplyr::select(-Entrez_Gene_Id)
all_dat <- all_dat[!(is.na(all_dat$Hugo_Symbol)),]
all_dat <- all_dat %>% group_by(Hugo_Symbol) %>% summarise_all(mean)
all_dat <- all_dat[match(inital_list, all_dat$Hugo_Symbol),]

for(i in 2:length(human_expre)){
  dat <- human_expre[i][[1]]
  dat <- dat[,-2]
  dat <- dat[!(is.na(dat$Hugo_Symbol)),]
  dat1 <- dat %>% group_by(Hugo_Symbol) %>% summarise_all(mean)
  dat1 <- dat1[match(inital_list, dat1$Hugo_Symbol),]
  message(paste0("BOUND ",i))
  print(dim(all_dat))
  print(all.equal(target = all_dat$Hugo_Symbol, current = dat1$Hugo_Symbol))
  if(all.equal(target = all_dat$Hugo_Symbol, current = dat1$Hugo_Symbol)){
    dat1 <- dat1 %>% dplyr::select(-Hugo_Symbol)
    all_dat <- cbind(all_dat, dat1)
  }
}
write.table(all_dat,here("data","homo_sapiens","panCancer_Atlas","final_all_panCancerExpression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

