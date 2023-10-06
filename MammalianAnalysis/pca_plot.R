## PCA plot for all dog samples ## 
# Load libraries
require(tibble) # dataframe manipulation
require(readr)  # for read_csv()
require(dplyr)  # for mutate()
require(tidyr)  # for unnest()
require(purrr)  # for map(), reduce()
library(here)

##########################################################INITIAL PCA ALL SAMPLES#########################################################################
## make PCA plot for only canine samples 
metadata <- readr::read_tsv(here("data","old_processed_data","canine_human_merged_metadata.tsv"))
canine_expr <- readr::read_tsv(here("data","processed_data","canine_logged_expression.tsv"))
canine_expr <- column_to_rownames(canine_expr, "Gene")
all_transposed <- t(canine_expr)
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
library(ggplot2)
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
canine_metadata <- metadata[metadata$Species == "Canine",]
canine_metadata[1:5,1:5]
dtp <- data.frame('Proj' = canine_metadata$Batch, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Proj)) + geom_point() + ggtitle("Canine Breed PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave("05_Canine_no_batch.png",width = 14, height = 8, dpi = 300)


## make PCA plot for only human samples 
human_expr <- readr::read_tsv(here("data","processed_data","human_logged_expression.tsv"))
human_expr <- column_to_rownames(human_expr, "Gene")
all_transposed <- t(human_expr)
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
human_metadata <- metadata[metadata$Species == "Human",]
human_metadata[1:5,1:5]
dtp <- data.frame('Proj' = human_metadata$Batch, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Proj)) + geom_point() + ggtitle("Human PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave("05_Human_no_batch.png",width = 14, height = 8, dpi = 300)

# both canine and human 
human_expr <- rownames_to_column(human_expr, "Gene")
canine_expr <- rownames_to_column(canine_expr, "Gene")
all_expression <- merge(human_expr, canine_expr, by = "Gene")
dim(all_expression)
all_expression[1:5,1:5]

#get rid of the normal samples from canine data 
metadata <- metadata[!(metadata$Tumor == "Normal"),]
all_expression <- column_to_rownames(all_expression, "Gene")
all_expression <- all_expression[,colnames(all_expression) %in% metadata$Sample]
all_transposed <- t(all_expression)
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Cancer' = metadata$CancerType, 'Species' = metadata$Species, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Proj)) + geom_point() + ggtitle("Human and Canine PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave("06_Human_Canine_no_batch_Cancer.png",width = 14, height = 8, dpi = 300)

dim(canine_metadata)
## here
# all_expression <- column_to_rownames(all, "Hugo_Symbol")
# all_logged <- log2(all_expression + 1)


##################BATCH CORRECTION FOR ALL SAMPLES ####################
#create batch ids
canine_batch <- readr::read_tsv(here("data","old_processed_data","canine_batch_ids.tsv"))
human_batch <- readr::read_tsv(here("data","old_processed_data","human_batch_ids.tsv"))

human <- human_batch$Batch
names(human) <- human_batch$Sample
canine <- canine_batch$Batch
names(canine) <- canine_batch$Sample
head(canine)
all_batch_ids <- c(human, canine)

#LESLIE TRY
#prepare for combat 
raw_merged <- as.matrix(all_expression) - min(all_expression)
dim(raw_merged)
library(sva)
dat_batch_adjusted_norm_new <- ComBat_seq(raw_merged, all_batch)
dat_batch_adjusted_norm_new[1:5,1:5]
#dat_batch_adjusted_norm_new <- as.data.frame(dat_batch_adjusted_norm_new)
#to_write <- rownames_to_column(dat_batch_adjusted_norm_new, "Gene")
#write.table(to_write, here("data","processed_data","batch_corrected_human_canine_expression.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

# create PCA 
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
# dat_pca$x[1:5,1:5]
# dat <- as.data.frame(dat_pca$x)
#to_write_pca <- rownames_to_column(dat, "Sample")
#write.table(to_write_pca, here("data","processed_data","batch_corrected_human_canine_PCA.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)
# dat <- readr::read_tsv(here("data","processed_data","batch_corrected_human_canine_PCA.tsv"))
# dat <- column_to_rownames(dat,"Sample")
#this section is messy currently from running things in termianl and rstudio - needs to be cleaned up:
#dtp <- data.frame("Cancer"=metadata$CancerType, 'Species' = metadata$Species, "dat" = dat[,1:2])
dtp <- data.frame("Cancer"=metadata$CancerType, 'Species' = metadata$Species, "dat" = dat_pca$x[,1:2])
#dtp <- data.frame("Cancer"=metadata$CancerType, 'Species' = metadata$Species, "dat" = all_pca$x[,1:2])
dtp$Species <- as.factor(dtp$Species)
levels(dtp$Species)

dat <- readr::read_tsv(here("data","processed_data","canine_human_PCA.tsv"))
dat <- column_to_rownames(dat,"Sample")
canine_metadata <- readr::read_tsv(here("data","canis_familiaris","metadata_final.txt"))
human_metadata <- readr::read_tsv(here("data","old_processed_data","final_human_metadata.txt"))

canine_cut <- canine_metadata %>% dplyr::select(Run, BioProject, Tumor, CancerType, Breed)
human_cut <- human_metadata %>% dplyr::select(Sample, ONCOTREE_CODE, TUMOR_TYPE,CANCER_TYPE)
colnames(canine_cut) <- c("Sample", "ProjCode","Tumor","CancerType", "Breed")
colnames(human_cut) <- c("Sample", "ProjCode","Tumor","CancerType")
#canine_cut$Breed <- "Canine"
human_cut$Breed<- "Human"

all_metadata <- rbind(human_cut, canine_cut)
all_metadata <- all_metadata[all_metadata$Tumor != "Normal",]
all_metadata <- all_metadata[all_metadata$Sample %in% rownames(dat),]
dim(all_metadata)

all_metadata$Breed[grep("Great Pyranees", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Great Pyrenees"


all_metadata$Breed[grep("Pitbull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
all_metadata$Breed[grep("PitBull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
all_metadata$Breed[grep("Labrador", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Labrador Retriever"
all_metadata$Breed[grep("Labrador mix", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Labrador Retriever"
all_metadata$Breed[grep("Labrador Retriever Mix", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Labrador Retriever"
all_metadata$Breed[grep("West Highland White terrier", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "West Highland White Terrier"

all_metadata$Breed[grep("Mix", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Mix Breed"
all_metadata$Breed[grep("mixed", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Mix Breed"
all_metadata$Breed[grep("King charles spaniel", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "King Charles Spaniel"
all_metadata$Breed[grep("French Bull Dog", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "French Bulldog"
all_metadata$Breed[grep("French bulldog", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "French Bulldog"
all_metadata$Breed[grep("Boston Terr", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Boston Terrier"
all_metadata$Breed[grep("Boston", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Boston Terrier"

all_metadata$Breed[grep("Yorkshire terrier", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Yorkshire Terrier"
all_metadata$Breed[grep("Scottish terrier", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Scottish Terrier"
all_metadata$Breed[grep("German Shepherd Dog", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "German Shepherd"
all_metadata$Breed[grep("Beagle X", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Beagle"
all_metadata$Breed[grep("Pomerian", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pomeranian"
all_metadata$Breed[grep("PitBull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
all_metadata$Breed[grep("PitBull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
all_metadata$Breed[grep("PitBull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
all_metadata$Breed[grep("PitBull", all_metadata$Breed, ignore.case = FALSE, fixed = TRUE)] <- "Pit Bull"
v

namet <- all_metadata[is.na(all_metadata$Breed),]

head(namet)







dtp <- data.frame("Cancer"=metadata$CancerType, 'Species' = metadata$Species, "dat" = dat[,1:2])
library(ggplot2)
q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color = Cancer)) + geom_point()+
  ggtitle("Human and Dog Cancer ") + labs(y= "PC2", x = "PC1") #+scale_color_manual(values = c("#b80904","#5A5A5A"))
q
ggsave("06human_canine_batch_corrected_PCA_CombatSeq_ColoredSpecies.png",width = 14, height = 8, dpi = 300)
# ###END 
# cale_shape_manual(values=c(3, 16, 17))
# scale_color_manual(values=c('#999999','#E69F00', '#56B4E9'))
# # - trasnpose matrix for PCA function - #
# dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)

# - get rid of 0 values - 0
# dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
# set.seed(417)
# dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
# dtp <- data.frame("Cancer"=all_metadata$CANCER_TYPE_DETAILED, 'Species' = all_metadata$SPECIES, "dat" = dat_pca$x[,1:2])
# q <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2, color = Cancer, shape=Species)) + geom_point() + ggtitle("Human and Dog Cancer Batch Correction") + labs(y= "PC2", x = "PC1")
# q
# #ggsave("MY_Test_Canine_human_batch_correction1.png",width = 30, height = 21, dpi = 300)
# ggsave("MY_Test_Canine_human_batch_correction1.png",width = 25, height = 16, dpi = 300)
# ggsave("MY_Test_Canine_human_batch_correction2.png",width = 14, height = 8, dpi = 300)
# ggsave("LeslieHelpCorrection.png",width = 14, height = 8, dpi = 300)

# dat_batch_adjusted_norm_new <- as.data.frame(dat_batch_adjusted_norm_new)
# to_write <- rownames_to_column(dat_batch_adjusted_norm_new, "Gene")
# write.table(to_write, here("data","batch_corrected_human_canine.tsv"),sep = "\t", col.names = TRUE, row.names = FALSE)






