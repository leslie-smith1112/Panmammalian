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
metadata <- readr::read_tsv(here("data","processed_data","canine_human_merged_metadata.tsv"))
canine_expr <- readr::read_tsv(here("data","processed_data","canine_expr.tsv"))
canine_expr <- column_to_rownames(canine_expr, "Gene")
all_transposed <- t(canine_expr)

# inital PCA for canine data 
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
library(ggplot2)
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
canine_metadata <- metadata[metadata$Species == "Canine",]
canine_metadata[1:5,1:5]
dtp <- data.frame('Proj' = canine_metadata$Batch, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Proj)) + geom_point() + ggtitle("Canine Breed PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave(here("plots","no_batch_bioproject.png"),width = 14, height = 8, dpi = 300)

#inial PCA for human data
human_expr <- readr::read_tsv(here("data","homo_sapiens","panCancer_Atlas","final_all_panCancerExpression.tsv"))
human_expr <- column_to_rownames(human_expr, "Gene")
all_transposed <- t(human_expr)
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
human_metadata <- metadata[metadata$Species == "Human",]
human_metadata[1:5,1:5]
dtp <- data.frame('Proj' = human_metadata$Batch, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Proj)) + geom_point() + ggtitle("Human PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave(here("plots","Human_Disease.png"),width = 14, height = 8, dpi = 300)

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

#initial canine and human PCA
all_transposed <- t(all_expression)
all_pcaMat <- all_transposed[ , which(apply(all_transposed, 2, var) != 0)]
all_pca <- prcomp(all_pcaMat, center = TRUE, scale. = TRUE)
dtp <- data.frame('Cancer' = metadata$CancerType, 'Species' = metadata$Species, dat = all_pca$x[,1:2])
p <- ggplot(data=dtp, aes(x=dat.PC1, y=dat.PC2,  color=Species)) + geom_point() + ggtitle("Human and Canine PCA - No Batch Correction") + labs(y= "PC2", x = "PC1")
p
ggsave("Humanan_Canina_Disease_PCA_no_Batch.png",width = 14, height = 8, dpi = 300)

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
all.equal(names(all_batch_ids), colnames(all_expression))
#prepare for combat 
raw_merged <- as.matrix(all_expression) - min(all_expression)
dim(raw_merged)
library(sva)
dat_batch_adjusted_norm_new <- ComBat_seq(raw_merged, all_batch)
dat_batch_adjusted_norm_new[1:5,1:5]
dat_batch_adjusted_norm_new <- as.data.frame(dat_batch_adjusted_norm_new)
to_write <- rownames_to_column(dat_batch_adjusted_norm_new, "Gene")
write.table(to_write,here("data","processed_data","batch_corrected_canine_human.tsv"), sep = "\t", col.names=TRUE, row.names = FALSE)

# create PCA 
dat_adjusted_norm_transposed <- t(dat_batch_adjusted_norm_new)
dat_adjusted_norm_pcaMat <- dat_adjusted_norm_transposed[ , which(apply(dat_adjusted_norm_transposed, 2, var) != 0)]
set.seed(417)
dat_pca <- prcomp(dat_adjusted_norm_pcaMat, center = TRUE, scale. = TRUE)
dat_pca$x[1:5,1:5]
dat <- as.data.frame(dat_pca$x)
to_write_pca <- rownames_to_column(dat, "Sample")
write.table(to_write_pca, here("data","processed_data","canine_human_PCA.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

dtp <- data.frame("Cancer"=metadata$CancerType, 'Species' = metadata$Species, "dat" = dat_pca$x[,1:2])
dtp$Species <- as.factor(dtp$Species)
levels(dtp$Species)

dat <- column_to_rownames(dat, "Sample")
metadata <- metadata[metadata$Sample %in% rownames(dat),]
library(ggplot2)
q <- ggplot(data=dtp, aes(x=dtp$dat.PC2, y=dtp$dat.PC3, color = Cancer)) + geom_point()+
  ggtitle("Human and Dog Cancer ") + labs(y= "PC2", x = "PC1") #+scale_color_manual(values = c("#b80904","#5A5A5A"))
q
ggsave("Canine_human_batch_corrected.png",width = 14, height = 8, dpi = 300)


#messing around with simplifying dog breeds - not used 
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
namet <- all_metadata[is.na(all_metadata$Breed),]







