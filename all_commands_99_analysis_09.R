###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	9.  KEGG overlap
#   - 9.1 overlap of genes in outlier windows between all four lakes
#   - 9.2 overlap of fst outlier windows with genes and the KEGG pathways of those genes
#   - 9.3 table output for paper

#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

library(tidyverse)

################# 9. 1 gene overlap from outlier windows ############################### 
#load genes that fall within outliers in each lake fst comparison
brienz_outlier_genes <- read.csv("KEGG/Brienz_extractedannotation_unique_gene_names.txt", header = FALSE)
neuchatel_outlier_genes <- read.csv("KEGG/Neuchatel_extractedannotation_unique_gene_names.txt", header = FALSE)
walen_outlier_genes <- read.csv("KEGG/Walen_extractedannotation_unique_gene_names.txt", header = FALSE)
luzern_outlier_genes <- read.csv("KEGG/Lucerne_extractedannotation_unique_gene_names.txt", header = FALSE)

brienz_outlier_genes$V1 <- as.character(brienz_outlier_genes$V1)
brienz_outlier_genes$lakeB <- "YES"
neuchatel_outlier_genes$V1 <- as.character(neuchatel_outlier_genes$V1)
neuchatel_outlier_genes$lakeN <- "YES"
walen_outlier_genes$V1 <- as.character(walen_outlier_genes$V1)
walen_outlier_genes$lakeW <- "YES"
luzern_outlier_genes$V1 <- as.character(luzern_outlier_genes$V1)
luzern_outlier_genes$lakeL <- "YES"

full_gene_df <- as.data.frame(full_join(brienz_outlier_genes, neuchatel_outlier_genes, by = "V1"))
full_gene_df <- as.data.frame(full_join(full_gene_df, walen_outlier_genes, by = "V1"))
full_gene_df <- as.data.frame(full_join(full_gene_df, luzern_outlier_genes, by = "V1"))

full_gene_df$count_empty <- rowSums(apply(is.na(full_gene_df), 2, as.numeric))
full_gene_df_shared <- subset(full_gene_df, full_gene_df$count_empty < 1)

################# 9. 2 KEGG overlap between comparisons ############################### 

#and load gene-kegg file
gene_to_KEGG <- read.csv("KEGG/all.maker.proteins.user_ko.txt", sep = "\t", header = FALSE)

#remove genes with no kegg annotation
gene_to_KEGG$empty <- "yes"
for (i in 1:nrow(gene_to_KEGG)){
  if(nchar(as.character(gene_to_KEGG$V2[i]))>0){
    gene_to_KEGG$empty[i] <- "no"
  }
}

#and make new df with only kegg annotated genes
gene_to_KEGG_clean <- subset(gene_to_KEGG, gene_to_KEGG$empty == "no")

gene_to_KEGG_clean$V1 <- as.character(gene_to_KEGG_clean$V1)
gene_to_KEGG_clean$V2 <- as.character(gene_to_KEGG_clean$V2)

gene_to_KEGG_clean$V1 <- gsub("-mRNA-1", "", gene_to_KEGG_clean$V1)

############First match the gene list we have with the gene lists for outliers in each lake###########
#make new empty vectors
brienz_vec <- c()
neuchatel_vec <- c()
walen_vec <- c()
luzern_vec <- c()

#and run loop so that if gene is in outlier window from that lake then it adds the kegg id of that gene into the column
for(i in 1:nrow(gene_to_KEGG_clean)){
  KEGG_ID_to_check <- as.character(gene_to_KEGG_clean$V1)[i]
  ID <- as.character(gene_to_KEGG_clean$V2)[i]
  a <- which(as.character(KEGG_ID_to_check) == as.vector(brienz_outlier_genes$V1))
  b <- which(as.character(KEGG_ID_to_check) == as.vector(neuchatel_outlier_genes$V1))
  c <- which(as.character(KEGG_ID_to_check) == as.vector(walen_outlier_genes$V1))
  d <- which(as.character(KEGG_ID_to_check) == as.vector(luzern_outlier_genes$V1))
  if(length(a)>0){
    brienz_vec[i] <- ID
  }
  if(length(a)<1){
    brienz_vec[i] <- "no"
  }
  if(length(b)>0){
    neuchatel_vec[i] <- ID
  }
  if(length(b)<1){
    neuchatel_vec[i] <- "no"
  }
  if(length(c)>0){
    walen_vec[i] <- ID
  }
  if(length(c)<1){
    walen_vec[i] <- "no"
  }
  if(length(d)>0){
    luzern_vec[i] <- ID
  }
  if(length(d)<1){
    luzern_vec[i] <- "no"
  }
}

gene_to_KEGG_clean$brienz <- brienz_vec
gene_to_KEGG_clean$neuchatel <- neuchatel_vec
gene_to_KEGG_clean$walen <- walen_vec
gene_to_KEGG_clean$luzern <- luzern_vec

#now sum number of missing values across rows
gene_to_KEGG_clean$count <- rowSums(gene_to_KEGG_clean[-1] == "no")

total_genes_present_at_least_once_across_4_lakesDF <- subset(gene_to_KEGG_clean, gene_to_KEGG_clean$count <5)
total_genes_present_at_least_once_across_4_lakes <- nrow(total_genes_present_at_least_once_across_4_lakesDF)
total_genes_present_all_4_lakesDF <- subset(gene_to_KEGG_clean, gene_to_KEGG_clean$count <2)
total_genes_present_all_4_lakes <- nrow(total_genes_present_all_4_lakesDF)
#so now we have four vectors with the KEGG IDs associated with outliers in each lake (some may occur more than once)

###########now use the list of KEGG IDs in each lake outlier to see if some are shared across all four###########
#now get unique kegg list from original df
unique_KEGGs <- unique(gene_to_KEGG_clean[c("V2")])

brienz_PRESENT <- c()
neuchatel_PRESENT <- c()
walen_PRESENT <- c()
luzern_PRESENT <- c()

for(i in 1:nrow(unique_KEGGs)){
  K_ID <- as.character(unique_KEGGs$V2)[i]
  aa <- which(as.character(K_ID) == brienz_vec)
  bb <- which(as.character(K_ID) == neuchatel_vec)
  cc <- which(as.character(K_ID) == walen_vec)
  dd <- which(as.character(K_ID) == luzern_vec)
  if(length(aa)>0){
    brienz_PRESENT[i] <- "yes"
  }
  if(length(aa)<1){
    brienz_PRESENT[i] <- "no"
  }
  if(length(bb)>0){
    neuchatel_PRESENT[i] <- "yes"
  }
  if(length(bb)<1){
    neuchatel_PRESENT[i] <- "no"
  }
  if(length(cc)>0){
    walen_PRESENT[i] <- "yes"
  }
  if(length(cc)<1){
    walen_PRESENT[i] <- "no"
  }
  if(length(dd)>0){
    luzern_PRESENT[i] <- "yes"
  }
  if(length(dd)<1){
    luzern_PRESENT[i] <- "no"
  }
}

FULL_KEGG_ID <- as.data.frame(unique_KEGGs)
FULL_KEGG_ID$brienz_ID <- brienz_PRESENT
FULL_KEGG_ID$neuchatel_ID <- neuchatel_PRESENT
FULL_KEGG_ID$walen_ID <- walen_PRESENT
FULL_KEGG_ID$luzern_ID <- luzern_PRESENT

#and count yes across columns
FULL_KEGG_ID$count <- rowSums(FULL_KEGG_ID[-1] == "yes")

colnames(FULL_KEGG_ID) <- c("KEGG_ID","brienz_ID","neuchatel_ID","walen_ID","luzern_ID","count")   
write.csv(FULL_KEGG_ID, "KEGG/FULL_KEGG_ID.csv", row.names=FALSE)

#calculate and write all kegg ids present at least once
morethan1 <- subset(FULL_KEGG_ID, FULL_KEGG_ID$count > 0)
write.csv(morethan1, "KEGG/morethan1.csv", row.names=FALSE)

KEGGID_SHARED_ALL4 <- subset(morethan1, morethan1$count > 3)

#write all lake-specific kegg id lists
brienz_KEGGs <- subset(morethan1, morethan1$brienz_ID == "yes")
neuchatel_KEGGs <- subset(morethan1, morethan1$neuchatel_ID == "yes")
walen_KEGGs <- subset(morethan1, morethan1$walen_ID == "yes")
luzern_KEGGs <- subset(morethan1, morethan1$luzern_ID == "yes")

write.csv(brienz_KEGGs, "KEGG/brienz_KEGGs.csv", row.names=FALSE)
write.csv(neuchatel_KEGGs, "KEGG/neuchatel_KEGGs.csv", row.names=FALSE)
write.csv(walen_KEGGs, "KEGG/walen_KEGGs.csv", row.names=FALSE)
write.csv(luzern_KEGGs, "KEGG/luzern_KEGGs.csv", row.names=FALSE)

#now the KEGG IDs are pasted into: https://www.kegg.jp/kegg/ko.html and the output copied back into an excel spreadsheet
#and read back in here:

B_KO <- read.csv("KEGG/brienz_KEGG_ID_OUTPUT.csv", header = FALSE, na.strings=c(""))
B_KO <- na.omit(B_KO)
colnames(B_KO) <- c("V1", "B_path", "B_count")

N_KO <- read.csv("KEGG/neuchatel_KEGG_ID_OUTPUT.csv", header = FALSE, na.strings=c(""))
N_KO <- N_KO[,1:3]
N_KO <- na.omit(N_KO)
colnames(N_KO) <- c("V1", "N_path", "N_count")

W_KO <- read.csv("KEGG/walen_KEGG_ID_OUTPUT.csv", header = FALSE, na.strings=c(""))
W_KO <- W_KO[,1:3]
W_KO <- na.omit(W_KO)
colnames(W_KO) <- c("V1", "W_path", "W_count")

L_KO <- read.csv("KEGG/luzern_KEGG_ID_OUTPUT.csv", header = FALSE, na.strings=c(""))
L_KO <- L_KO[,1:3]
L_KO <- na.omit(L_KO)
colnames(L_KO) <- c("V1", "L_path", "L_count")

nrow(B_KO)
nrow(N_KO)
nrow(W_KO)
nrow(L_KO)

#and join them all together to get a list of shared pathways
shared_KEGG_pathwaysDF <- full_join(B_KO, N_KO, by = "V1")
shared_KEGG_pathwaysDF <- full_join(shared_KEGG_pathwaysDF, W_KO, by = "V1")
shared_KEGG_pathwaysDF <- full_join(shared_KEGG_pathwaysDF, L_KO, by = "V1")
shared_KEGG_pathwaysDF$count_empty <- rowSums(apply(is.na(shared_KEGG_pathwaysDF), 2, as.numeric))

shared_KEGG_pathways <- subset(shared_KEGG_pathwaysDF, shared_KEGG_pathwaysDF$count_empty < 1)

################# 9. 3 table output for paper ############################### 
#table output:
##Genes overlapping with FST outlier windows
#total:
nrow(full_gene_df)
#shared:
nrow(full_gene_df_shared)

##Genes with associated KEGG IDs
#total = 
total_genes_present_at_least_once_across_4_lakes
#shared = 
total_genes_present_all_4_lakes

##Unique KEGG IDs from these genes
#total = 
nrow(morethan1)
#shared = 
nrow(KEGGID_SHARED_ALL4)

#KEGG pathways from these KEGG IDs
#total =
nrow(shared_KEGG_pathwaysDF)
#shared =
nrow(shared_KEGG_pathways)



