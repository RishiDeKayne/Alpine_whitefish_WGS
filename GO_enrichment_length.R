#response to reviewers
setwd("Dropbox/RishiMAC/99_reanalysis/")
#do outlive genes have the same length distribution as all other genes
GOraw <- read.csv("Parallel/GO_filt.out", sep = " ", stringsAsFactors = F)
outliers <- read.table("reviewer_responses/gene_names1659.txt",header=F,col.names = c("gene"))

all_genes_functional_names <- as.data.frame(unique(GOraw$qpid))
outliers_names <- outliers

all_gene_lengths <- read.csv("../../../Downloads/allr3.genes_filt_out.txt", sep = " ", head = F)
all_gene_lengths$length <- (all_gene_lengths$V3 - all_gene_lengths$V2)
all_gene_lengths$V1 <- as.character(all_gene_lengths$V1)
all_gene_lengths$V1 <- paste(all_gene_lengths$V1, "-mRNA-1", sep = "")


all_genes_functional_names$`unique(GOraw$qpid)` <- as.character(all_genes_functional_names$`unique(GOraw$qpid)`)
all_genes_functional_names$length <- 0
for (i in 1:nrow(all_genes_functional_names)){
  name <- all_genes_functional_names$`unique(GOraw$qpid)`[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    all_genes_functional_names$length[i] <- name_df$length
  }
}

outliers_names$gene <- as.character(outliers_names$gene)
outliers_names$gene <- paste(outliers_names$gene, "-mRNA-1", sep = "")

outliers_names$length <- 0
for (i in 1:nrow(outliers_names)){
  name <- outliers_names$gene[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    outliers_names$length[i] <- name_df$length
  }
}

all_genes_functional_names_len <- subset(all_genes_functional_names, all_genes_functional_names$length > 0)
outliers_names_len <- subset(outliers_names, outliers_names$length > 0)

all <- all_genes_functional_names_len$length
outliers <- outliers_names_len$length

par(mfrow=c(2,1))
hist(all)
hist(outliers)

par(mfrow=c(1,2))
boxplot(all, outliers, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("Genes overlapping CSS outlier windows n=", length(outliers), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)
wilcox.test(all,outliers)

library(vioplot)
#vioplot(all, outliers)

all_df <- as.data.frame(all)
all_df$numb <- 1
colnames(all_df) <- c("length", "numb")
outliers_df <- as.data.frame(outliers)
outliers_df$numb <- 2
colnames(outliers_df) <- c("length", "numb")
all_outliers_df <- rbind(all_df, outliers_df)


summary(lm(all_outliers_df$numb ~ all_outliers_df$length))

outliers1 <- outliers
########
#same for kegg
l_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Lucerne_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
l_genes$lake <- "L"
b_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Brienz_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
b_genes$lake <- "B"
w_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Walen_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
w_genes$lake <- "W"
n_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Neuchatel_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
n_genes$lake <- "N"

all_genes <- rbind(l_genes, b_genes, w_genes, n_genes)

all_genes$shared <- "no"
for (i in 1:nrow(all_genes)){
  gene_term_name <- as.character(all_genes$gene)[i]
  gene_term_name_df <- subset(all_genes, as.character(all_genes$gene) == gene_term_name)
  shared_by <- nrow(gene_term_name_df)
  all_genes$shared[i] <- shared_by
}

all_genes_unique <- all_genes[!duplicated(all_genes$gene),]
#write.table(all_genes_unique,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique.csv",sep=",",col.names = T,row.names = F)
all_genes_unique_shared <- subset(all_genes_unique, as.numeric(all_genes_unique$shared) > 1)
#write.table(all_genes_unique_shared,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique_shared.csv",sep=",",col.names = T,row.names = F)

kegg_genes <- as.data.frame(all_genes_unique$gene)
colnames(kegg_genes) <- c("gene")

outliers <- kegg_genes

all_genes_functional_names <- as.data.frame(unique(GOraw$qpid))
outliers_names <- outliers

all_gene_lengths <- read.csv("../../../Downloads/allr3.genes_filt_out.txt", sep = " ", head = F)
all_gene_lengths$length <- (all_gene_lengths$V3 - all_gene_lengths$V2)
all_gene_lengths$V1 <- as.character(all_gene_lengths$V1)
all_gene_lengths$V1 <- paste(all_gene_lengths$V1, "-mRNA-1", sep = "")


all_genes_functional_names$`unique(GOraw$qpid)` <- as.character(all_genes_functional_names$`unique(GOraw$qpid)`)
all_genes_functional_names$length <- 0
for (i in 1:nrow(all_genes_functional_names)){
  name <- all_genes_functional_names$`unique(GOraw$qpid)`[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    all_genes_functional_names$length[i] <- as.numeric(name_df$length)
  }
}

outliers_names$gene <- as.character(outliers_names$gene)
outliers_names$gene <- paste(outliers_names$gene, "-mRNA-1", sep = "")

outliers_names$length <- 0
for (i in 1:nrow(outliers_names)){
  name <- outliers_names$gene[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    outliers_names$length[i] <- name_df$length
  }
}

all_genes_functional_names_len <- subset(all_genes_functional_names, all_genes_functional_names$length > 0)
outliers_names_len <- subset(outliers_names, outliers_names$length > 0)

all <- all_genes_functional_names_len$length
outliers <- outliers_names_len$length

#par(mfrow=c(2,1))
#hist(all)
#hist(outliers)

boxplot(all, outliers, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("Genes overlapping FST outlier windows n=", length(outliers), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)
wilcox.test(all,outliers)

library(vioplot)
#vioplot(all, outliers)

all_df <- as.data.frame(all)
all_df$numb <- 1
colnames(all_df) <- c("length", "numb")
outliers_df <- as.data.frame(outliers)
outliers_df$numb <- 2
colnames(outliers_df) <- c("length", "numb")
all_outliers_df <- rbind(all_df, outliers_df)


summary(lm(all_outliers_df$numb ~ all_outliers_df$length))

outliers2 <- outliers

#plot
dev.off()
par(mfrow=c(1,2))
boxplot(all, outliers1, outliers2, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("CSS outlier genes n=", length(outliers1), sep = "")), (paste("FST outlier genes n=", length(outliers2), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)
vioplot(all, outliers1, outliers2, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("CSS outlier genes n=", length(outliers1), sep = "")), (paste("FST outlier genes n=", length(outliers2), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)


################################################
#download /cluster/work/gdc/shared/p659/analysis/99_reanalysis/parallel/Parallel_CSS_outliers/scaffoldannotation2.bed
#response to reviewers
setwd("Dropbox/RishiMAC/99_reanalysis/")
#do outlive genes have the same length distribution as all other genes
outliers <- read.table("reviewer_responses/gene_names1659.txt",header=F,col.names = c("gene"))

outliers <- unique(outliers)

all_gene_lengths <- read.csv("../../../Downloads/allr3.genes_filt_out.txt", sep = " ", head = F)
all_gene_lengths$length <- (all_gene_lengths$V3 - all_gene_lengths$V2)
all_gene_lengths$V1 <- as.character(all_gene_lengths$V1)
all_gene_lengths$V1 <- paste(all_gene_lengths$V1, "-mRNA-1", sep = "")

#the following file has all 42695 genes which were checked against outliers i.e AED filtered and on scaffolds
annotation <- read.table("reviewer_responses/scaffoldannotation2.bed",header=F)
annotation_split <- data.frame(do.call('rbind', strsplit(as.character(annotation$V4),';Name=',fixed=TRUE)))
annotation_names <- as.data.frame(annotation_split$X2)
colnames(annotation_names) <- c("gene")

annotation_split2 <- data.frame(do.call('rbind', strsplit(as.character(annotation_names$gene),';pseudo=',fixed=TRUE)))
annotation_names2 <- as.data.frame(annotation_split2$X1)
colnames(annotation_names2) <- c("gene")

annotation_names2$length <- 0
annotation_names2$gene <- paste(annotation_names2$gene, "-mRNA-1", sep = "")

for (i in 1:nrow(annotation_names2)){
  name <- annotation_names2$gene[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    annotation_names2$length[i] <- name_df$length
  }
}

outliers_names <- outliers
outliers_names$gene <- as.character(outliers_names$gene)
outliers_names$gene <- paste(outliers_names$gene, "-mRNA-1", sep = "")

outliers_names$length <- 0
for (i in 1:nrow(outliers_names)){
  name <- outliers_names$gene[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    outliers_names$length[i] <- name_df$length
  }
}


all <- annotation_names2$length
css_outliers <- outliers_names
outliers <- outliers_names$length

par(mfrow=c(2,1))
hist(all, breaks = 100)
hist(outliers, breaks = 50)

par(mfrow=c(1,2))
boxplot(all, outliers, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("Genes overlapping CSS outlier windows n=", length(outliers), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)


outliers1 <- outliers
########
#same for kegg
l_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Lucerne_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
l_genes$lake <- "L"
b_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Brienz_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
b_genes$lake <- "B"
w_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Walen_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
w_genes$lake <- "W"
n_genes <- read.table("Parallel/Parallel_fst/GO_outliers/Neuchatel_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))
n_genes$lake <- "N"

all_genes <- rbind(l_genes, b_genes, w_genes, n_genes)

all_genes$shared <- "no"
for (i in 1:nrow(all_genes)){
  gene_term_name <- as.character(all_genes$gene)[i]
  gene_term_name_df <- subset(all_genes, as.character(all_genes$gene) == gene_term_name)
  shared_by <- nrow(gene_term_name_df)
  all_genes$shared[i] <- shared_by
}

all_genes_unique <- all_genes[!duplicated(all_genes$gene),]
#write.table(all_genes_unique,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique.csv",sep=",",col.names = T,row.names = F)
all_genes_unique_shared <- subset(all_genes_unique, as.numeric(all_genes_unique$shared) > 1)
#write.table(all_genes_unique_shared,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique_shared.csv",sep=",",col.names = T,row.names = F)

kegg_genes <- as.data.frame(all_genes_unique$gene)
colnames(kegg_genes) <- c("gene")

outliers <- kegg_genes

outliers_names <- outliers

outliers_names$gene <- as.character(outliers_names$gene)
outliers_names$gene <- paste(outliers_names$gene, "-mRNA-1", sep = "")

outliers_names$length <- 0
for (i in 1:nrow(outliers_names)){
  name <- outliers_names$gene[i]
  name_df <- subset(all_gene_lengths, as.character(all_gene_lengths$V1) == name)
  if (nrow(name_df)>0){
    outliers_names$length[i] <- name_df$length
  }
}

outliers <- outliers_names$length
fst_outliers <- outliers_names


boxplot(all, outliers, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("Genes overlapping FST outlier windows n=", length(outliers), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5)

outliers2 <- outliers

#plot
dev.off()
par(mfrow=c(1,1))
boxplot(all, outliers1, outliers2, 
        names = c((paste("All genes n=", length(all), sep = "")), (paste("CSS outlier genes n=", length(outliers1), sep = "")), (paste("FST outlier genes n=", length(outliers2), sep = ""))),
        ylab = "gene length (bp)", cex = 0.5, col = "white")
#vioplot(all, outliers1, outliers2, 
#        names = c((paste("All genes n=", length(all), sep = "")), (paste("CSS outlier genes n=", length(outliers1), sep = "")), (paste("FST outlier genes n=", length(outliers2), sep = ""))),
#        ylab = "gene length (bp)", cex = 0.5)


#test statistical differences
#all = annotation_names2
#css = css_outliers
#fst = fst_outliers

#first test all vs. css
annotation_names2$ccs_overlap <- "no"
for(i in 1:nrow(annotation_names2)){
  gene_name <- as.character(annotation_names2$gene[i])
  css_test_df <- subset(css_outliers, as.character(css_outliers$gene) == gene_name)
  if(nrow(css_test_df)>0){
    annotation_names2$ccs_overlap[i] <- "yes"
  }
}

annotation_withcssyes <- subset(annotation_names2, annotation_names2$ccs_overlap == "yes")
annotation_withcssno <- subset(annotation_names2, annotation_names2$ccs_overlap == "no")

annotation_names2$fst_overlap <- "no"
for(i in 1:nrow(annotation_names2)){
  gene_name <- as.character(annotation_names2$gene[i])
  fst_test_df <- subset(fst_outliers, as.character(fst_outliers$gene) == gene_name)
  if(nrow(fst_test_df)>0){
    annotation_names2$fst_overlap[i] <- "yes"
  }
}

annotation_withcssyes <- subset(annotation_names2, annotation_names2$ccs_overlap == "yes")
annotation_withcssno <- subset(annotation_names2, annotation_names2$ccs_overlap == "no")

annotation_withfstyes <- subset(annotation_names2, annotation_names2$fst_overlap == "yes")
annotation_withfstno <- subset(annotation_names2, annotation_names2$fst_overlap == "no")

wilcox.test(annotation_withcssno$length, outliers1, alternative = "less")
#W = 32178234, p-value = 2.775e-08
nrow(annotation_withcssno) #40993
length(outliers1) #1702

wilcox.test(annotation_withfstno$length, outliers2, alternative = "less")
#W = 21624555, p-value = 2.693e-06
nrow(annotation_withfstno) #41565
length(outliers2) #1130
boxplot(annotation_withcssno$length, outliers1, annotation_withfstno$length, outliers2,
        ylab = "gene length (bp)", cex = 0.5, col = "white")

hist(all)
hist(outliers1)
hist(outliers2)

p1 <- hist(all)                     # centered at 4
p2 <- hist(outliers1)                     # centered at 6
p3 <- hist(outliers2)                     # centered at 6

plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,150000))  # first histogra m
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,150000), add=T)
plot( p3, col=rgb(0,1,0,1/4), xlim=c(0,150000), add=T)

plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,150000))
plot( p3, col=rgb(0,1,0,1/4), xlim=c(0,150000), add=T)

