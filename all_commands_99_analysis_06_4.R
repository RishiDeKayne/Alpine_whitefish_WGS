###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	6.  parallel
#see other script for   - 6.1 parallel PCA
#see other script for		- 6.2 Fst
#see other script for	- 6.3 CSS
#		- 6.4 outlier analysis
#     -6.4.1 GO enrichment of 342 windows
#     -6.4.2 PCA with 342 windows
#     -6.4.3 PCA with 342 windows with outgroups
#     -6.4.4 PCA with 342 windows without gill raker associated chromosome scaff22

setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

################# 6.4.1 GO enrichment of 342 windows ############################### 
# GO analysis
library(topGO)

# set node size
nodeS <- 10

# .tsv with 6 columns: "locus.name","relationship","GO.term","GO.ID","aspect","evidence.code"
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
GOraw <- read.csv("Parallel/GO_filt.out", sep = " ", stringsAsFactors = F)
#outliers <- read.table("./gene_names.txt",header=F,col.names = c("gene"))
outliers <- read.table("Parallel/gene_names342.txt",header=F,col.names = c("gene"))

outliers$gene <- paste(outliers$gene, "-mRNA-1", sep = "")
outliers$gene <- as.factor(outliers$gene)

# split into the three ontologies
GOraw_C <- subset(GOraw, ontology == "CC")
GOraw_B <- subset(GOraw, ontology == "BP")
GOraw_M <- subset(GOraw, ontology == "MF")


# make gene2GO objects for each ontology (a list where each element is a vector of all GO.IDs for a gebe, named with the gene name)
mygene2GO_C_raw <- sapply(GOraw_C$qpid,function(x) as.character(unique(GOraw_C$GO.goid[GOraw_C$qpid==x])))
mygene2GO_B_raw <- sapply(GOraw_B$qpid,function(x) as.character(unique(GOraw_B$GO.goid[GOraw_B$qpid==x])))
mygene2GO_M_raw <- sapply(GOraw_M$qpid,function(x) as.character(unique(GOraw_M$GO.goid[GOraw_M$qpid==x])))

GOraw_C$gene_factor <- as.factor(GOraw_C$qpid)
GOraw_B$gene_factor <- as.factor(GOraw_B$qpid)
GOraw_M$gene_factor <- as.factor(GOraw_M$qpid)

names(mygene2GO_C_raw) <- GOraw_C$gene_factor
names(mygene2GO_B_raw) <- GOraw_B$gene_factor
names(mygene2GO_M_raw) <- GOraw_M$gene_factor


# make geneList (factor of all genes with 1=outlier and 0=not) 
geneNames <- unique(GOraw_C$gene_factor)
sigGenes <- outliers$gene
geneList <- factor(as.integer(geneNames %in% sigGenes))
names(geneList) <- geneNames

# make TopGO data objects for each ontology
GOdata_C_raw <- new("topGOdata",
                    ontology="CC",
                    allGenes=geneList,
                    annot=annFUN.gene2GO,
                    gene2GO=mygene2GO_C_raw,
                    nodeSize = nodeS)

GOdata_B_raw <- new("topGOdata",
                    ontology="BP",
                    allGenes=geneList,
                    annot=annFUN.gene2GO,
                    gene2GO=mygene2GO_B_raw,
                    nodeSize = nodeS)

GOdata_M_raw <- new("topGOdata",
                    ontology="MF",
                    allGenes=geneList,
                    annot=annFUN.gene2GO,
                    gene2GO=mygene2GO_M_raw,
                    nodeSize = nodeS)

# Run tests
Fisher.weig.CC_raw <- runTest(GOdata_C_raw, algorithm = "weight", statistic = "fisher")
Fisher.elim.CC_raw <- runTest(GOdata_C_raw, algorithm = "elim", statistic = "fisher")
Fisher.clas.CC_raw <- runTest(GOdata_C_raw, algorithm = "classic", statistic = "fisher")
res.CC_raw = GenTable(GOdata_C_raw, weight_fisher_P = Fisher.weig.CC_raw, elim_fisher_P = Fisher.elim.CC_raw, classic_fisher_P = Fisher.clas.CC_raw, topNodes = 375, numChar = 1000)
res.CC_raw_uncorrected <- subset(res.CC_raw, Significant > 1 & (as.numeric(weight_fisher_P) < 0.05 & as.numeric(elim_fisher_P) < 0.05))

allGO_CC = usedGO(object = GOdata_C_raw) 
# use it in GenTable as follows:
res.CC_cor = GenTable(GOdata_C_raw, weight_fisher_P = Fisher.weig.CC_raw, elim_fisher_P = Fisher.elim.CC_raw, classic_fisher_P = Fisher.clas.CC_raw, topNodes = length(allGO_CC), numChar = 1000)

res.CC_cor$weight_fisher_P_cor <- p.adjust(res.CC_cor$weight_fisher_P, method = "fdr", n = length(res.CC_cor$weight_fisher_P))
res.CC_cor$elim_fisher_P_cor <- p.adjust(res.CC_cor$elim_fisher_P, method = "fdr", n = length(res.CC_cor$elim_fisher_P))
res.CC_cor$classic_fisher_P_cor <- p.adjust(res.CC_cor$classic_fisher_P, method = "fdr", n = length(res.CC_cor$classic_fisher_P))
res.CC_cor_fdrcorrected <- subset(res.CC_cor, Significant > 1 & (as.numeric(weight_fisher_P_cor) < 0.05 & as.numeric(elim_fisher_P_cor) < 0.05))


Fisher.weig.BP_raw <- runTest(GOdata_B_raw, algorithm = "weight", statistic = "fisher")
Fisher.elim.BP_raw <- runTest(GOdata_B_raw, algorithm = "elim", statistic = "fisher")
Fisher.clas.BP_raw <- runTest(GOdata_B_raw, algorithm = "classic", statistic = "fisher")
res.BP_raw = GenTable(GOdata_B_raw, weight_fisher_P = Fisher.weig.BP_raw, elim_fisher_P = Fisher.elim.BP_raw, classic_fisher_P = Fisher.clas.BP_raw, topNodes = 375, numChar = 1000)
res.BP_raw_uncorrected <- subset(res.BP_raw, Significant > 1 & (as.numeric(weight_fisher_P) < 0.05 & as.numeric(elim_fisher_P) < 0.05))

allGO_BP = usedGO(object = GOdata_B_raw) 
# use it in GenTable as follows:
res.BP_cor = GenTable(GOdata_B_raw, weight_fisher_P = Fisher.weig.BP_raw, elim_fisher_P = Fisher.elim.BP_raw, classic_fisher_P = Fisher.clas.BP_raw, topNodes = length(allGO_BP), numChar = 1000)

res.BP_cor$weight_fisher_P_cor <- p.adjust(res.BP_cor$weight_fisher_P, method = "fdr", n = length(res.BP_cor$weight_fisher_P))
res.BP_cor$elim_fisher_P_cor <- p.adjust(res.BP_cor$elim_fisher_P, method = "fdr", n = length(res.BP_cor$elim_fisher_P))
res.BP_cor$classic_fisher_P_cor <- p.adjust(res.BP_cor$classic_fisher_P, method = "fdr", n = length(res.BP_cor$classic_fisher_P))
res.BP_caw_fdrcorrected <- subset(res.BP_cor, Significant > 1 & (as.numeric(weight_fisher_P_cor) < 0.05 & as.numeric(elim_fisher_P_cor) < 0.05))


Fisher.weig.MF_raw <- runTest(GOdata_M_raw, algorithm = "weight", statistic = "fisher")
Fisher.elim.MF_raw <- runTest(GOdata_M_raw, algorithm = "elim", statistic = "fisher")
Fisher.clas.MF_raw <- runTest(GOdata_M_raw, algorithm = "classic", statistic = "fisher")
res.MF_raw = GenTable(GOdata_M_raw, weight_fisher_P = Fisher.weig.MF_raw, elim_fisher_P = Fisher.elim.MF_raw, classic_fisher_P = Fisher.clas.MF_raw, topNodes = 375, numChar = 1000)
res.MF_raw_uncorrected <- subset(res.MF_raw, Significant > 1 & (as.numeric(weight_fisher_P) < 0.05 & as.numeric(elim_fisher_P) < 0.05))

allGO_MF = usedGO(object = GOdata_M_raw) 
# use it in GenTable as follows:
res.MF_cor = GenTable(GOdata_M_raw, weight_fisher_P = Fisher.weig.MF_raw, elim_fisher_P = Fisher.elim.MF_raw, classic_fisher_P = Fisher.clas.MF_raw, topNodes = length(allGO_MF), numChar = 1000)

res.MF_cor$weight_fisher_P_cor <- p.adjust(res.MF_cor$weight_fisher_P, method = "fdr", n = length(res.MF_cor$weight_fisher_P))
res.MF_cor$elim_fisher_P_cor <- p.adjust(res.MF_cor$elim_fisher_P, method = "fdr", n = length(res.MF_cor$elim_fisher_P))
res.MF_cor$classic_fisher_P_cor <- p.adjust(res.MF_cor$classic_fisher_P, method = "fdr", n = length(res.MF_cor$classic_fisher_P))
res.MF_cor_fdrcorrected <- subset(res.MF_cor, Significant > 1 & (as.numeric(weight_fisher_P_cor) < 0.05 & as.numeric(elim_fisher_P_cor) < 0.05))


# print results
write.table(res.CC_raw,file=paste0("Parallel/Parallel_css/GOenrichment_CC_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.BP_raw,file=paste0("Parallel/Parallel_css/GOenrichment_BP_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.MF_raw,file=paste0("Parallel/Parallel_css/GOenrichment_MF_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)

CC_output <- res.CC_raw_uncorrected
CC_output$class <- "CC"

BP_output <- res.BP_raw_uncorrected
BP_output$class <- "BP"

MF_output <- res.MF_raw_uncorrected
MF_output$class <- "MF"

full_GO_output <- rbind(CC_output, BP_output, MF_output)
write.table(full_GO_output,file="Parallel/Parallel_css/full_GO_output_nodes10.csv",sep=",",col.names = T,row.names = F)

################# 6.4.2 CSS outlier windows PCA ############################### 

# sort out the pca data
##11526 variants and 91 people pass filters and QC.

#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca <- read_table2("Parallel/All_342_allsnps_filt.eigenvec", col_names = FALSE)
eigenval <- scan("Parallel/All_342_allsnps_filt.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
for (i in 1:nrow(pca)){
  x <- pca$X2[i]
  pca$X2[i] <- substring(x, 2)
}

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))

pca$ind <- as.integer(pca$ind)

#get species from background
for (i in 1:nrow(pca)){
  speciesname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  spp[i] <- as.character(speciesname$Spp_mod)
}

# location
loc <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  locname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  loc[i] <- as.character(locname$Lake)
}

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

fishec_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  fishec_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  fishec_plot[i] <- as.character(fishec_name$fishec_ID)
}

#add colour
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$lake_col)
}

#add ecomorph
eco_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  eco_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  eco_plot[i] <- as.character(eco_name$ecomorph_mod)
}

# remake data.frame
#pca <- as.tibble(data.frame(pca, spp, loc, spp_loc, fishec_plot))
pca <- as.data.frame(data.frame(pca, spp, loc, spp_loc, fishec_plot, colour_plot, eco_plot))

#add symbol
pca$symb <- c()
for (i in 1:nrow(pca)){
  if(pca$eco_plot[i] == "balchen"){
    pca$symb[i] <- 22
  }
  if(pca$eco_plot[i] == "albeli"){
    pca$symb[i] <- 21
  }
  if(pca$eco_plot[i] == "felchen"){
    pca$symb[i] <- 23
  }
  if(pca$eco_plot[i] == "large_pelagic"){
    pca$symb[i] <- 8
  }
  if(pca$eco_plot[i] == "pelagic_profundal"){
    pca$symb[i] <- 24
  }
  if(pca$eco_plot[i] == "benthic_profundal"){
    pca$symb[i] <- 25
  }
}

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

par(mfrow=c(1,1))

#
#tiff("Parallel/Parallel_css/All_342PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.05)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P2")
#text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
#dev.off()

plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P3")
#text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))

plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC2)-0.005), (max(pca$PC2)+0.1)),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P2/P3")
#text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))

Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lty = 2)
} 


hull_balchen <- subset(pca, pca$eco_plot == "balchen")
hull_albeli <- subset(pca, pca$eco_plot == "albeli")
hull_felchen <- subset(pca, pca$eco_plot == "felchen")
hull_LP <- subset(pca, pca$eco_plot == "large_pelagic")
hull_PP <- subset(pca, pca$eco_plot == "pelagic_profundal")
hull_BP <- subset(pca, pca$eco_plot == "benthic_profundal")

#tiff("Parallel/Parallel_css/All_342PC1PC2_hulls.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("Parallel/Parallel_css/All_342PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)), (max(pca$PC1))),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P2")
#dev.off()
Plot_ConvexHull(xcoord = hull_balchen$PC1, ycoord = hull_balchen$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_albeli$PC1, ycoord = hull_albeli$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_felchen$PC1, ycoord = hull_felchen$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_LP$PC1, ycoord = hull_LP$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_PP$PC1, ycoord = hull_PP$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_BP$PC1, ycoord = hull_BP$PC2, lcolor = "black")
#dev.off()

newplot_df <- subset(pca, pca$loc == "Thun" | pca$loc == "Lucerne")
newplot_df$plotnumber <- c()

for (i in 1:nrow(newplot_df)){
  if (newplot_df$loc[i] == "Thun"){
    newplot_df$plotnumber[i] <- 0.1}
  if (newplot_df$loc[i] == "Lucerne"){
    newplot_df$plotnumber[i] <- 0.2}
  if (newplot_df$eco_plot[i] == "large_pelagic" | newplot_df$eco_plot[i] == "pelagic_profundal" | newplot_df$eco_plot[i] == "benthic_profundal"){
    newplot_df$plotnumber[i] <- newplot_df$plotnumber[i]-0.05} 
}

par(mar=c(5,9,4,1)+.1)
plot(newplot_df$PC1, newplot_df$plotnumber,
     col = as.character(newplot_df$colour_plot),
     pch = newplot_df$symb,
     lwd = 2,
     ylim = c(0,0.3),
     cex = 2, 
     yaxt='n',
     ylab= "",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
axis(2, at = c(0.05, 0.1, 0.15, 0.2), labels=c("Thun - Rare", "Thun - Common", "Lucerne - Rare", "Lucerne - Common"), las=1, cex = 1)

newplot2_df <- pca
newplot2_df$plotnumber <- c()

for (i in 1:nrow(newplot2_df)){
  if (newplot2_df$loc[i] == "Constance"){
    newplot2_df$plotnumber[i] <- 0.1}
  if (newplot2_df$loc[i] == "Lucerne"){
    newplot2_df$plotnumber[i] <- 0.2}
  if (newplot2_df$loc[i] == "Thun"){
    newplot2_df$plotnumber[i] <- 0.3}
  if (newplot2_df$eco_plot[i] == "large_pelagic" | newplot2_df$eco_plot[i] == "pelagic_profundal" | newplot2_df$eco_plot[i] == "benthic_profundal"){
    newplot2_df$plotnumber[i] <- newplot2_df$plotnumber[i]-0.05} 
  if (newplot2_df$loc[i] == "Brienz"){
    newplot2_df$plotnumber[i] <- 0.35}
  if (newplot2_df$loc[i] == "Zurich"){
    newplot2_df$plotnumber[i] <- 0.4}
  if (newplot2_df$loc[i] == "Walen"){
    newplot2_df$plotnumber[i] <- 0.45}
  if (newplot2_df$loc[i] == "Biel"){
    newplot2_df$plotnumber[i] <- 0.5}
  if (newplot2_df$loc[i] == "Neuenburg"){
    newplot2_df$plotnumber[i] <- 0.55}
}

#plot lakewise PC1 vs lake
#tiff("Parallel/Parallel_css/All_342PC1PC2_by_lake.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("Parallel/Parallel_css/All_342PC1PC2_by_lake_hulls.tiff", height=8, width=8, units="in", res=300, compression="lzw")
tiff("Parallel/Parallel_css/All_342PC1_vs_lake.tiff", height=8, width=12, units="in", res=300, compression="lzw")
par(mar=c(5,9,4,1)+.2)
plot(newplot2_df$PC1, newplot2_df$plotnumber,
     col = as.character(newplot2_df$colour_plot),
     pch = newplot2_df$symb,
     lwd = 2,
     ylim = c(0,0.6),
     cex = 2, 
     yaxt='n',
     ylab= "",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
axis(2, at = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55), labels=c("Constance - Rare", "Constance - Common", "Lucerne - Rare", "Lucerne - Common", "Thun - Rare", "Thun - Common", "Brienz", "Zurich", "Walen", "Biel", "Neuenburg"), las=1, cex = 1)
abline(h=0.025, col = "black", lty = 3)
abline(h=0.125, col = "black", lty = 3)
abline(h=0.225, col = "black", lty = 3)
abline(h=0.375, col = "black", lty = 3)
abline(h=0.475, col = "black", lty = 3)
abline(h=0.575, col = "black", lty = 3)
#text(newplot2_df$PC1, newplot2_df$plotnumber, newplot2_df$ind, cex = 0.4, pos = 3)
dev.off()

hull_balchen <- subset(newplot2_df, newplot2_df$eco_plot == "balchen")
hull_albeli <- subset(newplot2_df, newplot2_df$eco_plot == "albeli")
hull_felchen <- subset(newplot2_df, newplot2_df$eco_plot == "felchen")
hull_LP <- subset(newplot2_df, newplot2_df$eco_plot == "large_pelagic")
hull_PP <- subset(newplot2_df, newplot2_df$eco_plot == "pelagic_profundal")
hull_BP <- subset(newplot2_df, newplot2_df$eco_plot == "benthic_profundal")

Plot_ConvexHull(xcoord = hull_balchen$PC1, ycoord = hull_balchen$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_albeli$PC1, ycoord = hull_albeli$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_felchen$PC1, ycoord = hull_felchen$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_LP$PC1, ycoord = hull_LP$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_PP$PC1, ycoord = hull_PP$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_BP$PC1, ycoord = hull_BP$plotnumber, lcolor = "black")
#dev.off()


#########################################################################
#plots to look at correlation of gill raker/Standard length with PC1
par(mar=c(5,4,4,1)+.1)

gill <- pca

gill$gill_plot <- c()
for (i in 1:nrow(gill)){
  indiv_name <- as.character(gill$ind)[i]
  gill_count_sub <- subset(background, as.character(background$Indiv) == indiv_name)
  gill$gill_plot[i] <- as.character(gill_count_sub$gill_raker_count)
}

gill <- subset(gill, gill$gill_plot != "missing")
gill$gill_plot <- as.numeric(gill$gill_plot)

PC1_GRC <- lm(gill$gill_plot ~ gill$PC1)
summary(PC1_GRC)
#Multiple R-squared:  0.3734,	Adjusted R-squared:  0.3663 
#F-statistic: 52.45 on 1 and 88 DF,  p-value: 1.589e-10

GRC_r2 <- summary(PC1_GRC)$r.squared
GRC_r2_2dp <- format(round(GRC_r2, 3), nsmall = 3)


plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. gill raker count - R2 =", GRC_r2_2dp))
#abline(lm(gill$gill_plot ~ gill$PC1), lty = 3)

#and new the same but without C. profundus:
gill_noprof <- subset(gill, gill$spp != "C.profundus")
gill_noprof$gill_plot <- as.numeric(gill_noprof$gill_plot)

PC1_GRC_NOPROF <- lm(gill_noprof$gill_plot ~ gill_noprof$PC1)
summary(PC1_GRC_NOPROF)
#Multiple R-squared:  0.4337,	Adjusted R-squared:  0.427 
#F-statistic:  65.1 on 1 and 85 DF,  p-value: 4.124e-12

GRC_r2_NOPROF <- summary(PC1_GRC_NOPROF)$r.squared
GRC_r2_2dp_NOPROF <- format(round(GRC_r2_NOPROF, 3), nsmall = 3)

plot(gill_noprof$PC1, gill_noprof$gill_plot, col = as.character(gill_noprof$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. gill raker count - R2 =", GRC_r2_2dp_NOPROF))
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 3)




#split up lm by lake

####Luzern
luzern_gill <- subset(gill, gill$loc == "Lucerne")
par(mar=c(5,4,4,1)+.1)

PC1_GRC_L <- lm(luzern_gill$gill_plot ~ luzern_gill$PC1)
summary(PC1_GRC_L)
#Multiple R-squared:  0.8113,	Adjusted R-squared:  0.7995 
#F-statistic: 68.79 on 1 and 16 DF,  p-value: 3.462e-07

GRC_L_r2 <- summary(PC1_GRC_L)$r.squared
GRC_L_r2_2dp <- format(round(GRC_L_r2, 3), nsmall = 3)


plot(luzern_gill$PC1, luzern_gill$gill_plot, col = as.character(luzern_gill$colour_plot), 
     pch=luzern_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("Luzern PC1 vs. gill raker count - R2 =", GRC_L_r2_2dp))
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")

####Walen/Zurich
walen_gill <- subset(gill, gill$loc == "Walen" | gill$loc == "Zurich")
PC1_GRC_W <- lm(walen_gill$gill_plot ~ walen_gill$PC1)
summary(PC1_GRC_W)
#Multiple R-squared:  0.7195,	Adjusted R-squared:  0.702 
#F-statistic: 41.04 on 1 and 16 DF,  p-value: 8.693e-06

GRC_W_r2 <- summary(PC1_GRC_W)$r.squared
GRC_W_r2_2dp <- format(round(GRC_W_r2, 3), nsmall = 3)


plot(walen_gill$PC1, walen_gill$gill_plot, col = as.character(walen_gill$colour_plot), 
     pch=walen_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("walen PC1 vs. gill raker count - R2 =", GRC_W_r2_2dp))
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")


####Brienz/Thun - without profundal
brienz_gill <- subset(gill, gill$loc == "Brienz" | gill$loc == "Thun")
brienz_gill <- subset(brienz_gill, as.character(brienz_gill$spp) != "C.profundus")
PC1_GRC_B <- lm(brienz_gill$gill_plot ~ brienz_gill$PC1)
summary(PC1_GRC_B)
#Multiple R-squared:  0.753,	Adjusted R-squared:  0.7442 
#F-statistic: 85.37 on 1 and 28 DF,  p-value: 5.357e-10

GRC_B_r2 <- summary(PC1_GRC_B)$r.squared
GRC_B_r2_2dp <- format(round(GRC_B_r2, 3), nsmall = 3)


plot(brienz_gill$PC1, brienz_gill$gill_plot, col = as.character(brienz_gill$colour_plot), 
     pch=brienz_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("brienz/thun PC1 vs. gill raker count - R2 =", GRC_B_r2_2dp))
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")

####Neuchatel/biel
neuchatel_gill <- subset(gill, gill$loc == "Neuenburg" | gill$loc == "Biel")
PC1_GRC_N <- lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1)
summary(PC1_GRC_N)
#Multiple R-squared:  0.6813,	Adjusted R-squared:  0.6459 
#F-statistic: 19.24 on 1 and 9 DF,  p-value: 0.001755

GRC_N_r2 <- summary(PC1_GRC_N)$r.squared
GRC_N_r2_2dp <- format(round(GRC_N_r2, 3), nsmall = 3)


plot(neuchatel_gill$PC1, neuchatel_gill$gill_plot, col = as.character(neuchatel_gill$colour_plot), 
     pch=neuchatel_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("neuchatel/thun PC1 vs. gill raker count - R2 =", GRC_N_r2_2dp))
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")


####constance
constance_gill <- subset(gill, gill$loc == "Constance")
PC1_GRC_C <- lm(constance_gill$gill_plot ~ constance_gill$PC1)
summary(PC1_GRC_C)
#Multiple R-squared:  0.4375,	Adjusted R-squared:  0.3671 
#F-statistic: 6.221 on 1 and 8 DF,  p-value: 0.03727

GRC_C_r2 <- summary(PC1_GRC_C)$r.squared
GRC_C_r2_2dp <- format(round(GRC_C_r2, 3), nsmall = 3)


plot(constance_gill$PC1, constance_gill$gill_plot, col = as.character(constance_gill$colour_plot), 
     pch=constance_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("constance/thun PC1 vs. gill raker count - R2 =", GRC_C_r2_2dp))
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")


#now plot full plot again with all ab lines separately
tiff("Parallel/Parallel_css/All_342PC1_vs_GillRaker_noprof.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. gill raker count"))
abline(lm(gill$gill_plot ~ gill$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 2, col = "grey", lwd = 2)
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")
legend('topright', 
       legend = c((paste("Overall - R2 =", GRC_r2_2dp)),
                  (paste("Overall (exl. profundus) - R2 =", GRC_r2_2dp_NOPROF)),
                  (paste("Brienz/Thun (exl. profundus) - R2 =", GRC_B_r2_2dp)), 
                  (paste("Lucerne - R2 =", GRC_L_r2_2dp)), 
                  (paste("Constance - R2 =", GRC_C_r2_2dp)), 
                  (paste("Neuchatel/Biel - R2 =", GRC_N_r2_2dp)), 
                  (paste("Walen/Zurich - R2 =", GRC_W_r2_2dp))),
       col = c("black",
               "grey",
               "chartreuse4", 
               "#ff4489", 
               "goldenrod", 
               "turquoise4", 
               "purple4"),
       pt.cex = 1.2, cex = 0.8, pch = "-")
dev.off()

##FOR FIG5
tiff("Parallel/Parallel_css/All_342PC1_vs_GillRaker_noprof_FIG1.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. gill raker count"))
abline(lm(gill$gill_plot ~ gill$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 2, col = "grey", lwd = 2)
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")
#legend('topright', 
#       legend = c((paste("Overall - R2 =", GRC_r2_2dp)),
#                  (paste("Overall (exl. profundus) - R2 =", GRC_r2_2dp_NOPROF)),
#                  (paste("Brienz/Thun (exl. profundus) - R2 =", GRC_B_r2_2dp)), 
#                  (paste("Lucerne - R2 =", GRC_L_r2_2dp)), 
#                  (paste("Constance - R2 =", GRC_C_r2_2dp)), 
#                  (paste("Neuchatel/Biel - R2 =", GRC_N_r2_2dp)), 
#                  (paste("Walen/Zurich - R2 =", GRC_B_r2_2dp))),
#       col = c("black",
#               "grey",
#               "chartreuse4", 
#               "#ff4489", 
#               "goldenrod", 
#               "turquoise4", 
#               "purple4"),
#       pt.cex = 1.2, cex = 0.8, pch = "-")
dev.off()

##################################################
#standard length similar 
par(mar=c(5,4,4,1)+.1)

SL <- pca

SL$SL_plot <- c()
for (i in 1:nrow(SL)){
  indiv_name <- as.character(SL$ind)[i]
  SL_sub <- subset(background, as.character(background$Indiv) == indiv_name)
  SL$SL_plot[i] <- as.character(SL_sub$standard_length)
}

SL <- subset(SL, SL$SL_plot != "missing")
SL$SL_plot <- as.numeric(SL$SL_plot)

PC1_SL <- lm(SL$SL_plot ~ SL$PC1)
summary(PC1_SL)
#Multiple R-squared:  0.3846,	Adjusted R-squared:  0.3777 
#F-statistic: 55.63 on 1 and 89 DF,  p-value: 5.527e-11

summary(PC1_SL)$r.squared

SL_r2 <- summary(PC1_SL)$r.squared
SL_r2_2dp <- format(round(SL_r2, 3), nsmall = 3)

plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "standard length",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. standard length - R2 =", SL_r2_2dp))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 3)

####Luzern
luzern_SL <- subset(SL, SL$loc == "Lucerne")
par(mar=c(5,4,4,1)+.1)

PC1_SL_L <- lm(luzern_SL$SL_plot ~ luzern_SL$PC1)
summary(PC1_SL_L)
#Multiple R-squared:  0.4546,	Adjusted R-squared:  0.4206 
#F-statistic: 13.34 on 1 and 16 DF,  p-value: 0.002149

SL_L_r2 <- summary(PC1_SL_L)$r.squared
SL_L_r2_2dp <- format(round(SL_L_r2, 3), nsmall = 3)


plot(luzern_SL$PC1, luzern_SL$SL_plot, col = as.character(luzern_SL$colour_plot), 
     pch=luzern_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("Luzern PC1 vs. SL - R2 =", SL_L_r2_2dp))
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")

####Walen/Zurich
walen_SL <- subset(SL, SL$loc == "Walen" | SL$loc == "Zurich")
PC1_SL_W <- lm(walen_SL$SL_plot ~ walen_SL$PC1)
summary(PC1_SL_W)
#Multiple R-squared:  0.6986,	Adjusted R-squared:  0.6798 
#F-statistic: 37.08 on 1 and 16 DF,  p-value: 1.564e-05

SL_W_r2 <- summary(PC1_SL_W)$r.squared
SL_W_r2_2dp <- format(round(SL_W_r2, 3), nsmall = 3)


plot(walen_SL$PC1, walen_SL$SL_plot, col = as.character(walen_SL$colour_plot), 
     pch=walen_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("walen PC1 vs. SL - R2 =", SL_W_r2_2dp))
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")


####Brienz/Thun - without profundal
brienz_SL <- subset(SL, SL$loc == "Brienz" | SL$loc == "Thun")
PC1_SL_B <- lm(brienz_SL$SL_plot ~ brienz_SL$PC1)
summary(PC1_SL_B)
#Multiple R-squared:  0.3913,	Adjusted R-squared:  0.3717 
#F-statistic: 19.93 on 1 and 31 DF,  p-value: 9.903e-05

SL_B_r2 <- summary(PC1_SL_B)$r.squared
SL_B_r2_2dp <- format(round(SL_B_r2, 3), nsmall = 3)


plot(brienz_SL$PC1, brienz_SL$SL_plot, col = as.character(brienz_SL$colour_plot), 
     pch=brienz_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("brienz/thun PC1 vs. SL - R2 =", SL_B_r2_2dp))
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")

####Neuchatel/biel
neuchatel_SL <- subset(SL, SL$loc == "Neuenburg" | SL$loc == "Biel")
PC1_SL_N <- lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1)
summary(PC1_SL_N)
#Multiple R-squared:  0.6771,	Adjusted R-squared:  0.6448 
#F-statistic: 20.97 on 1 and 10 DF,  p-value: 0.001012

SL_N_r2 <- summary(PC1_SL_N)$r.squared
SL_N_r2_2dp <- format(round(SL_N_r2, 3), nsmall = 3)


plot(neuchatel_SL$PC1, neuchatel_SL$SL_plot, col = as.character(neuchatel_SL$colour_plot), 
     pch=neuchatel_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("neuchatel/thun PC1 vs. SL - R2 =", SL_N_r2_2dp))
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")


####constance
constance_SL <- subset(SL, SL$loc == "Constance")
PC1_SL_C <- lm(constance_SL$SL_plot ~ constance_SL$PC1)
summary(PC1_SL_C)
#Multiple R-squared:  0.0618,	Adjusted R-squared:  -0.05547 
#F-statistic: 0.527 on 1 and 8 DF,  p-value: 0.4886

SL_C_r2 <- summary(PC1_SL_C)$r.squared
SL_C_r2_2dp <- format(round(SL_C_r2, 3), nsmall = 3)


plot(constance_SL$PC1, constance_SL$SL_plot, col = as.character(constance_SL$colour_plot), 
     pch=constance_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("constance/thun PC1 vs. SL - R2 =", SL_C_r2_2dp))
abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")


#now plot full plot again with all ab lines separately
#tiff("Parallel/Parallel_css/All_342PC1_vs_StandardLength.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. SL"))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")
#abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")
legend('topleft', 
       legend = c((paste("Overall - R2 =", SL_r2_2dp)),
                  (paste("Brienz/Thun - R2 =", SL_B_r2_2dp)), 
                  (paste("Lucerne - R2 =", SL_L_r2_2dp)), 
                  (paste("Constance (non-significant) - R2 =", SL_C_r2_2dp)), 
                  (paste("Neuchatel/Biel - R2 =", SL_N_r2_2dp)), 
                  (paste("Walen/Zurich - R2 =", SL_W_r2_2dp))),
       col = c("black",
               "chartreuse4", 
               "#ff4489", 
               "goldenrod", 
               "turquoise4", 
               "purple4"),
       pt.cex = 1.2, cex = 0.8, pch = "-")
#dev.off()

#FOR FIG5
#tiff("Parallel/Parallel_css/All_342PC1_vs_StandardLength.tiff", height=8, width=10, units="in", res=300, compression="lzw")
tiff("Parallel/Parallel_css/All_342PC1_vs_StandardLength_Fig5.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. SL"))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")
#abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")
#legend('topleft', 
#       legend = c((paste("Overall - R2 =", SL_r2_2dp)),
#                  (paste("Brienz/Thun - R2 =", SL_B_r2_2dp)), 
#                  (paste("Lucerne - R2 =", SL_L_r2_2dp)), 
#                  (paste("Neuchatel/Biel - R2 =", SL_N_r2_2dp)), 
#                  (paste("Constance (non-significant) - R2 =", SL_C_r2_2dp)), 
#                  (paste("Walen/Zurich - R2 =", SL_W_r2_2dp))),
#       col = c("black",
#               "chartreuse4", 
#               "#ff4489", 
#               "goldenrod", 
#               "turquoise4", 
#               "purple4"),
#       pt.cex = 1.2, cex = 0.8, pch = "-")
dev.off()

#add figure from 05.R
#now analyse the 7 outlier snps and how the allele frequencies differ for each ecomorph
balchen_df <- read.csv("GWAS/balchen_freq.frq", header = F, sep = "\t")
balchen_df <- balchen_df[-1,]
albeli_df <- read.csv("GWAS/albeli_freq.frq", header = F, sep = "\t")
albeli_df <- albeli_df[-1,]
felchen_df <- read.csv("GWAS/felchen_freq.frq", header = F, sep = "\t")
felchen_df <- felchen_df[-1,]
LP_df <- read.csv("GWAS/largepel_freq.frq", header = F, sep = "\t")
LP_df <- LP_df[-1,]
benth_prof_df <- read.csv("GWAS/benth_prof_freq.frq", header = F, sep = "\t")
benth_prof_df <-  benth_prof_df[-1,]
pel_prof_df <- read.csv("GWAS/pel_prof_freq.frq", header = F, sep = "\t")
pel_prof_df <- pel_prof_df[-1,]

indiv_with_GR <- subset(background, background$gill_raker_count != "missing")

balchen_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "balchen")
balchen_gr_mean <- mean(as.numeric(as.character(balchen_gr_bg$gill_raker_count)))
albeli_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "albeli")
albeli_gr_mean <- mean(as.numeric(as.character(albeli_gr_bg$gill_raker_count)))
felchen_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "felchen")
felchen_gr_mean <- mean(as.numeric(as.character(felchen_gr_bg$gill_raker_count)))
LP_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "large_pelagic")
LP_gr_mean <- mean(as.numeric(as.character(LP_gr_bg$gill_raker_count)))
BP_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "benthic_profundal")
BP_gr_mean <- mean(as.numeric(as.character(BP_gr_bg$gill_raker_count)))
PP_gr_bg <- subset(indiv_with_GR, indiv_with_GR$ecomorph_mod == "pelagic_profundal")
PP_gr_mean <- mean(as.numeric(as.character(PP_gr_bg$gill_raker_count)))

#find out gill raker count means
balchen_gr_mean
#2
albeli_gr_mean
#5
felchen_gr_mean
#3
LP_gr_mean
#4
BP_gr_mean
#1
PP_gr_mean
#6

balchen_df$order_val <- 2
balchen_df$ecomorph <- "balchen"

albeli_df$order_val <- 5
albeli_df$ecomorph <- "albeli"

felchen_df$order_val <- 3
felchen_df$ecomorph <- "felchen"

LP_df$order_val <- 4
LP_df$ecomorph <- "LP"

benth_prof_df$order_val <-  1
benth_prof_df$ecomorph <- "BP"

pel_prof_df$order_val <- 6
pel_prof_df$ecomorph <- "PP"

full_snp_comparison <- rbind(benth_prof_df, balchen_df, felchen_df, LP_df, albeli_df, pel_prof_df)

full_snp_comparison$V5 <- gsub("A:", "", full_snp_comparison$V5)
full_snp_comparison$V5 <- gsub("T:", "", full_snp_comparison$V5)
full_snp_comparison$V5 <- gsub("G:", "", full_snp_comparison$V5)
full_snp_comparison$V5 <- gsub("C:", "", full_snp_comparison$V5)

full_snp_comparison$V6 <- gsub("A:", "", full_snp_comparison$V6)
full_snp_comparison$V6 <- gsub("T:", "", full_snp_comparison$V6)
full_snp_comparison$V6 <- gsub("G:", "", full_snp_comparison$V6)
full_snp_comparison$V6 <- gsub("C:", "", full_snp_comparison$V6)

##plot figure with allele frequencies for the 7 outlier snps by ecomorph
par(mar=c(4,5,4,5)+0.1)
plot(full_snp_comparison$order_val, full_snp_comparison$V5, ylab = "reference allele frequency", xlab = "ecomorph", lwd = 2, xaxt='n', pch=c(25, 22, 23, 8, 21, 24), col = "black", ylim = c(0,1), cex = 1)
axis(1, at = c(1,2,3,4,5,6), labels=c("Benthic profundal", "Balchen", "Felchen", "Large pelagic", "Albeli", "Pelagic profundal"), las=1, cex = 1, cex.axis=0.7)
lines(full_snp_comparison$order_val,full_snp_comparison$V5, col = "grey50", lty = 3)
#and add the gill raker points
par(new = T)
plot(c(1,2,3,4,5,6), c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean), axes=F, xlab=NA, ylab=NA, cex=1, col = "orangered4", ylim = c(15, 40), pch=c(25, 22, 23, 8, 21, 24), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'mean gill raker count', col = "orangered4", cex = 1)

#and add the gill raker points
par(new = T)
plot(c(1,2,3,4,5,6), c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean), axes=F, xlab=NA, ylab=NA, cex=1, col = "orangered4", ylim = c(15, 40), pch=c(25, 22, 23, 8, 21, 24), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'mean gill raker count', col = "orangered4", cex = 1)

#plot for fig5
par(mfrow=c(2,2))
tiff("../FIGURE5.tiff", height=12, width=12, units="in", res=300, compression="lzw")
#tiff("FIGURE5rev.tiff", height=12, width=12, units="in", res=300, compression="lzw")

#plot1
par(mar=c(4,9.5,4,1)+.1)
plot(newplot2_df$PC1, newplot2_df$plotnumber,
     col = as.character(newplot2_df$colour_plot),
     pch = newplot2_df$symb,
     lwd = 2,
     ylim = c(0,0.6),
     cex = 2, 
     yaxt='n',
     ylab= "",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
axis(2, at = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55), labels=c("Constance - Rare", "Constance - Common", "Lucerne - Rare", "Lucerne - Common", "Thun - Rare", "Thun - Common", "Brienz", "Zurich", "Walen", "Biel", "Neuenburg"), las=1, cex = 1)
abline(h=0.025, col = "black", lty = 3)
abline(h=0.125, col = "black", lty = 3)
abline(h=0.225, col = "black", lty = 3)
abline(h=0.375, col = "black", lty = 3)
abline(h=0.475, col = "black", lty = 3)
abline(h=0.575, col = "black", lty = 3)

#plot2
par(mar=c(4,5,4,5)+0.1)
plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "Standard length",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
#main = ("PC1 vs. SL"))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")

#plot3
par(mar=c(4,9.5,4,1)+.1)
plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "Gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
     #main = ("PC1 vs. gill raker count")
abline(lm(gill$gill_plot ~ gill$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 2, col = "grey", lwd = 2)
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")

#plot4
par(mar=c(4,5,4,5)+0.1)
plot(full_snp_comparison$order_val, full_snp_comparison$V5, ylab = "reference allele frequency", xlab = "ecomorph", lwd = 2, xaxt='n', pch=c(25, 22, 23, 8, 21, 24), col = "black", ylim = c(0,1), cex = 1.5)
axis(1, at = c(1,2,3,4,5,6), labels=c("Benthic profundal", "Balchen", "Felchen", "Large pelagic", "Albeli", "Pelagic profundal"), las=1, cex = 2.5, cex.axis=1)
lines(full_snp_comparison$order_val,full_snp_comparison$V5, col = "grey50", lty = 3)
#and add the gill raker points
par(new = T)
plot(c(1,2,3,4,5,6), c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean), axes=F, xlab=NA, ylab=NA, cex=1.5, col = "orangered4", ylim = c(15, 40), pch=c(25, 22, 23, 8, 21, 24), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'mean gill raker count', col = "orangered4", cex = 1)

dev.off()
##################################################

################# 6.4.3 CSS outlier windows PCA with outgroups ############################### 
#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca <- read_table2("Parallel/All_342_allsnps_wout_filt.eigenvec", col_names = FALSE)
eigenval <- scan("Parallel/All_342_allsnps_wout_filt.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
for (i in 1:nrow(pca)){
  x <- pca$X2[i]
  pca$X2[i] <- substring(x, 2)
}

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))

pca$ind <- as.integer(pca$ind)

#get species from background
for (i in 1:nrow(pca)){
  speciesname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  spp[i] <- as.character(speciesname$Spp_mod)
}

# location
loc <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  locname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  loc[i] <- as.character(locname$Lake)
}

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

fishec_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  fishec_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  fishec_plot[i] <- as.character(fishec_name$fishec_ID)
}

#add colour
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$lake_col)
}

#add ecomorph
eco_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  eco_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  eco_plot[i] <- as.character(eco_name$ecomorph_mod)
}

# remake data.frame
#pca <- as.tibble(data.frame(pca, spp, loc, spp_loc, fishec_plot))
pca <- as.data.frame(data.frame(pca, spp, loc, spp_loc, fishec_plot, colour_plot, eco_plot))

#add symbol
pca$symb <- c()
for (i in 1:nrow(pca)){
  if(pca$eco_plot[i] == "balchen"){
    pca$symb[i] <- 22
  }
  if(pca$eco_plot[i] == "albeli"){
    pca$symb[i] <- 21
  }
  if(pca$eco_plot[i] == "felchen"){
    pca$symb[i] <- 23
  }
  if(pca$eco_plot[i] == "large_pelagic"){
    pca$symb[i] <- 8
  }
  if(pca$eco_plot[i] == "pelagic_profundal"){
    pca$symb[i] <- 24
  }
  if(pca$eco_plot[i] == "benthic_profundal"){
    pca$symb[i] <- 25
  }
}

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

par(mfrow=c(1,1))

#
#tiff("Parallel/Parallel_css/All_342PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.05)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples PCA with 342 CS' outlier windows with outgroups - PC1/PC2")
#text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
#dev.off()

plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P3")
#text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))

plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC2)-0.005), (max(pca$PC2)+0.1)),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P2/P3")
#text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))


################# 6.4.4 CSS outlier windows PCA without Gill raker count association peak ############################### 
#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca <- read_table2("Parallel/All_342_allsnps_noGRCpeak_filt.eigenvec", col_names = FALSE)
eigenval <- scan("Parallel/All_342_allsnps_noGRCpeak_filt.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
for (i in 1:nrow(pca)){
  x <- pca$X2[i]
  pca$X2[i] <- substring(x, 2)
}

# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))

pca$ind <- as.integer(pca$ind)

#get species from background
for (i in 1:nrow(pca)){
  speciesname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  spp[i] <- as.character(speciesname$Spp_mod)
}

# location
loc <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  locname <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  loc[i] <- as.character(locname$Lake)
}

# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

fishec_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  fishec_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  fishec_plot[i] <- as.character(fishec_name$fishec_ID)
}

#add colour
colour_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  col_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  colour_plot[i] <- as.character(col_name$lake_col)
}

#add ecomorph
eco_plot <- rep(NA, length(pca$ind))
for (i in 1:nrow(pca)){
  eco_name <- subset(background, as.character(background$Indiv) == as.character(pca$ind[i]))
  eco_plot[i] <- as.character(eco_name$ecomorph_mod)
}

# remake data.frame
#pca <- as.tibble(data.frame(pca, spp, loc, spp_loc, fishec_plot))
pca <- as.data.frame(data.frame(pca, spp, loc, spp_loc, fishec_plot, colour_plot, eco_plot))

#add symbol
pca$symb <- c()
for (i in 1:nrow(pca)){
  if(pca$eco_plot[i] == "balchen"){
    pca$symb[i] <- 22
  }
  if(pca$eco_plot[i] == "albeli"){
    pca$symb[i] <- 21
  }
  if(pca$eco_plot[i] == "felchen"){
    pca$symb[i] <- 23
  }
  if(pca$eco_plot[i] == "large_pelagic"){
    pca$symb[i] <- 8
  }
  if(pca$eco_plot[i] == "pelagic_profundal"){
    pca$symb[i] <- 24
  }
  if(pca$eco_plot[i] == "benthic_profundal"){
    pca$symb[i] <- 25
  }
}

# first convert to percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

par(mfrow=c(1,1))

#
#tiff("Parallel/Parallel_css/All_342PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.05)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - 342 CS' outlier windows - no scaff22 GRC peak")
#text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
#dev.off()

plot(pca$PC1, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P3")
#text(pca$PC1, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))

plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC2)-0.005), (max(pca$PC2)+0.1)),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P2/P3")
#text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))

Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  lines(xcoord[hpts], ycoord[hpts], col = lcolor, lty = 2)
} 


hull_balchen <- subset(pca, pca$eco_plot == "balchen")
hull_albeli <- subset(pca, pca$eco_plot == "albeli")
hull_felchen <- subset(pca, pca$eco_plot == "felchen")
hull_LP <- subset(pca, pca$eco_plot == "large_pelagic")
hull_PP <- subset(pca, pca$eco_plot == "pelagic_profundal")
hull_BP <- subset(pca, pca$eco_plot == "benthic_profundal")

#tiff("Parallel/Parallel_css/All_342PC1PC2_hulls.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("Parallel/Parallel_css/All_342PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)), (max(pca$PC1))),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P2")
#dev.off()
Plot_ConvexHull(xcoord = hull_balchen$PC1, ycoord = hull_balchen$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_albeli$PC1, ycoord = hull_albeli$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_felchen$PC1, ycoord = hull_felchen$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_LP$PC1, ycoord = hull_LP$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_PP$PC1, ycoord = hull_PP$PC2, lcolor = "black")
Plot_ConvexHull(xcoord = hull_BP$PC1, ycoord = hull_BP$PC2, lcolor = "black")
#dev.off()

newplot_df <- subset(pca, pca$loc == "Thun" | pca$loc == "Lucerne")
newplot_df$plotnumber <- c()

for (i in 1:nrow(newplot_df)){
  if (newplot_df$loc[i] == "Thun"){
    newplot_df$plotnumber[i] <- 0.1}
  if (newplot_df$loc[i] == "Lucerne"){
    newplot_df$plotnumber[i] <- 0.2}
  if (newplot_df$eco_plot[i] == "large_pelagic" | newplot_df$eco_plot[i] == "pelagic_profundal" | newplot_df$eco_plot[i] == "benthic_profundal"){
    newplot_df$plotnumber[i] <- newplot_df$plotnumber[i]-0.05} 
}

par(mar=c(5,9,4,1)+.1)
plot(newplot_df$PC1, newplot_df$plotnumber,
     col = as.character(newplot_df$colour_plot),
     pch = newplot_df$symb,
     lwd = 2,
     ylim = c(0,0.3),
     cex = 2, 
     yaxt='n',
     ylab= "",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")))
axis(2, at = c(0.05, 0.1, 0.15, 0.2), labels=c("Thun - Rare", "Thun - Common", "Lucerne - Rare", "Lucerne - Common"), las=1, cex = 1)

newplot2_df <- pca
newplot2_df$plotnumber <- c()

for (i in 1:nrow(newplot2_df)){
  if (newplot2_df$loc[i] == "Constance"){
    newplot2_df$plotnumber[i] <- 0.1}
  if (newplot2_df$loc[i] == "Lucerne"){
    newplot2_df$plotnumber[i] <- 0.2}
  if (newplot2_df$loc[i] == "Thun"){
    newplot2_df$plotnumber[i] <- 0.3}
  if (newplot2_df$eco_plot[i] == "large_pelagic" | newplot2_df$eco_plot[i] == "pelagic_profundal" | newplot2_df$eco_plot[i] == "benthic_profundal"){
    newplot2_df$plotnumber[i] <- newplot2_df$plotnumber[i]-0.05} 
  if (newplot2_df$loc[i] == "Brienz"){
    newplot2_df$plotnumber[i] <- 0.35}
  if (newplot2_df$loc[i] == "Zurich"){
    newplot2_df$plotnumber[i] <- 0.4}
  if (newplot2_df$loc[i] == "Walen"){
    newplot2_df$plotnumber[i] <- 0.45}
  if (newplot2_df$loc[i] == "Biel"){
    newplot2_df$plotnumber[i] <- 0.5}
  if (newplot2_df$loc[i] == "Neuenburg"){
    newplot2_df$plotnumber[i] <- 0.55}
}

#plot lakewise PC1 vs lake
#tiff("Parallel/Parallel_css/All_342PC1PC2_by_lake.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("Parallel/Parallel_css/All_342PC1PC2_by_lake_hulls.tiff", height=8, width=8, units="in", res=300, compression="lzw")
par(mar=c(5,9,4,1)+.1)
plot(newplot2_df$PC1, newplot2_df$plotnumber,
     col = as.character(newplot2_df$colour_plot),
     pch = newplot2_df$symb,
     lwd = 2,
     ylim = c(0,0.6),
     cex = 2, 
     yaxt='n',
     ylab= "",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = "CS' outlier windows no scaff22")
axis(2, at = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55), labels=c("Constance - Rare", "Constance - Common", "Lucerne - Rare", "Lucerne - Common", "Thun - Rare", "Thun - Common", "Brienz", "Zurich", "Walen", "Biel", "Neuenburg"), las=1, cex = 1)
abline(h=0.025, col = "black", lty = 3)
abline(h=0.125, col = "black", lty = 3)
abline(h=0.225, col = "black", lty = 3)
abline(h=0.375, col = "black", lty = 3)
abline(h=0.475, col = "black", lty = 3)
abline(h=0.575, col = "black", lty = 3)
#text(newplot2_df$PC1, newplot2_df$plotnumber, newplot2_df$ind, cex = 0.4, pos = 3)
#dev.off()

hull_balchen <- subset(newplot2_df, newplot2_df$eco_plot == "balchen")
hull_albeli <- subset(newplot2_df, newplot2_df$eco_plot == "albeli")
hull_felchen <- subset(newplot2_df, newplot2_df$eco_plot == "felchen")
hull_LP <- subset(newplot2_df, newplot2_df$eco_plot == "large_pelagic")
hull_PP <- subset(newplot2_df, newplot2_df$eco_plot == "pelagic_profundal")
hull_BP <- subset(newplot2_df, newplot2_df$eco_plot == "benthic_profundal")

Plot_ConvexHull(xcoord = hull_balchen$PC1, ycoord = hull_balchen$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_albeli$PC1, ycoord = hull_albeli$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_felchen$PC1, ycoord = hull_felchen$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_LP$PC1, ycoord = hull_LP$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_PP$PC1, ycoord = hull_PP$plotnumber, lcolor = "black")
Plot_ConvexHull(xcoord = hull_BP$PC1, ycoord = hull_BP$plotnumber, lcolor = "black")
#dev.off()


#########################################################################
#plots to look at correlation of gill raker/Standard length with PC1
par(mar=c(5,4,4,1)+.1)

gill <- pca

gill$gill_plot <- c()
for (i in 1:nrow(gill)){
  indiv_name <- as.character(gill$ind)[i]
  gill_count_sub <- subset(background, as.character(background$Indiv) == indiv_name)
  gill$gill_plot[i] <- as.character(gill_count_sub$gill_raker_count)
}

gill <- subset(gill, gill$gill_plot != "missing")
gill$gill_plot <- as.numeric(gill$gill_plot)

PC1_GRC <- lm(gill$gill_plot ~ gill$PC1)
summary(PC1_GRC)
#Multiple R-squared:  0.3734,	Adjusted R-squared:  0.3663 
#F-statistic: 52.45 on 1 and 88 DF,  p-value: 1.589e-10

GRC_r2 <- summary(PC1_GRC)$r.squared
GRC_r2_2dp <- format(round(GRC_r2, 3), nsmall = 3)


plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. gill raker count - R2 =", GRC_r2_2dp))
#abline(lm(gill$gill_plot ~ gill$PC1), lty = 3)

#and new the same but without C. profundus:
gill_noprof <- subset(gill, gill$spp != "C.profundus")
gill_noprof$gill_plot <- as.numeric(gill_noprof$gill_plot)

PC1_GRC_NOPROF <- lm(gill_noprof$gill_plot ~ gill_noprof$PC1)
summary(PC1_GRC_NOPROF)
#Multiple R-squared:  0.4337,	Adjusted R-squared:  0.427 
#F-statistic:  65.1 on 1 and 85 DF,  p-value: 4.124e-12

GRC_r2_NOPROF <- summary(PC1_GRC_NOPROF)$r.squared
GRC_r2_2dp_NOPROF <- format(round(GRC_r2_NOPROF, 3), nsmall = 3)

plot(gill_noprof$PC1, gill_noprof$gill_plot, col = as.character(gill_noprof$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. gill raker count - R2 =", GRC_r2_2dp_NOPROF))
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 3)




#split up lm by lake

####Luzern
luzern_gill <- subset(gill, gill$loc == "Lucerne")
par(mar=c(5,4,4,1)+.1)

PC1_GRC_L <- lm(luzern_gill$gill_plot ~ luzern_gill$PC1)
summary(PC1_GRC_L)
#Multiple R-squared:  0.8113,	Adjusted R-squared:  0.7995 
#F-statistic: 68.79 on 1 and 16 DF,  p-value: 3.462e-07

GRC_L_r2 <- summary(PC1_GRC_L)$r.squared
GRC_L_r2_2dp <- format(round(GRC_L_r2, 3), nsmall = 3)


plot(luzern_gill$PC1, luzern_gill$gill_plot, col = as.character(luzern_gill$colour_plot), 
     pch=luzern_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("Luzern PC1 vs. gill raker count - R2 =", GRC_L_r2_2dp))
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")

####Walen/Zurich
walen_gill <- subset(gill, gill$loc == "Walen" | gill$loc == "Zurich")
PC1_GRC_W <- lm(walen_gill$gill_plot ~ walen_gill$PC1)
summary(PC1_GRC_W)
#Multiple R-squared:  0.7195,	Adjusted R-squared:  0.702 
#F-statistic: 41.04 on 1 and 16 DF,  p-value: 8.693e-06

GRC_W_r2 <- summary(PC1_GRC_W)$r.squared
GRC_W_r2_2dp <- format(round(GRC_W_r2, 3), nsmall = 3)


plot(walen_gill$PC1, walen_gill$gill_plot, col = as.character(walen_gill$colour_plot), 
     pch=walen_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("walen PC1 vs. gill raker count - R2 =", GRC_W_r2_2dp))
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")


####Brienz/Thun - without profundal
brienz_gill <- subset(gill, gill$loc == "Brienz" | gill$loc == "Thun")
brienz_gill <- subset(brienz_gill, as.character(brienz_gill$spp) != "C.profundus")
PC1_GRC_B <- lm(brienz_gill$gill_plot ~ brienz_gill$PC1)
summary(PC1_GRC_B)
#Multiple R-squared:  0.753,	Adjusted R-squared:  0.7442 
#F-statistic: 85.37 on 1 and 28 DF,  p-value: 5.357e-10

GRC_B_r2 <- summary(PC1_GRC_B)$r.squared
GRC_B_r2_2dp <- format(round(GRC_B_r2, 3), nsmall = 3)


plot(brienz_gill$PC1, brienz_gill$gill_plot, col = as.character(brienz_gill$colour_plot), 
     pch=brienz_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("brienz/thun PC1 vs. gill raker count - R2 =", GRC_B_r2_2dp))
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")

####Neuchatel/biel
neuchatel_gill <- subset(gill, gill$loc == "Neuenburg" | gill$loc == "Biel")
PC1_GRC_N <- lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1)
summary(PC1_GRC_N)
#Multiple R-squared:  0.6813,	Adjusted R-squared:  0.6459 
#F-statistic: 19.24 on 1 and 9 DF,  p-value: 0.001755

GRC_N_r2 <- summary(PC1_GRC_N)$r.squared
GRC_N_r2_2dp <- format(round(GRC_N_r2, 3), nsmall = 3)


plot(neuchatel_gill$PC1, neuchatel_gill$gill_plot, col = as.character(neuchatel_gill$colour_plot), 
     pch=neuchatel_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("neuchatel/thun PC1 vs. gill raker count - R2 =", GRC_N_r2_2dp))
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")


####constance
constance_gill <- subset(gill, gill$loc == "Constance")
PC1_GRC_C <- lm(constance_gill$gill_plot ~ constance_gill$PC1)
summary(PC1_GRC_C)
#Multiple R-squared:  0.4375,	Adjusted R-squared:  0.3671 
#F-statistic: 6.221 on 1 and 8 DF,  p-value: 0.03727

GRC_C_r2 <- summary(PC1_GRC_C)$r.squared
GRC_C_r2_2dp <- format(round(GRC_C_r2, 3), nsmall = 3)


plot(constance_gill$PC1, constance_gill$gill_plot, col = as.character(constance_gill$colour_plot), 
     pch=constance_gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("constance/thun PC1 vs. gill raker count - R2 =", GRC_C_r2_2dp))
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")


#now plot full plot again with all ab lines separately
#tiff("Parallel/Parallel_css/All_342PC1_vs_GillRaker_noscaff22_noprof.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. gill raker count"))
abline(lm(gill$gill_plot ~ gill$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(gill_noprof$gill_plot ~ gill_noprof$PC1), lty = 2, col = "grey", lwd = 2)
abline(lm(luzern_gill$gill_plot ~ luzern_gill$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_gill$gill_plot ~ walen_gill$PC1), lty = 3, col = "purple4")
abline(lm(brienz_gill$gill_plot ~ brienz_gill$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_gill$gill_plot ~ neuchatel_gill$PC1), lty = 3, col = "turquoise4")
abline(lm(constance_gill$gill_plot ~ constance_gill$PC1), lty = 3, col = "goldenrod")
legend('topright', 
       legend = c((paste("Overall - R2 =", GRC_r2_2dp)),
                  (paste("Overall (exl. profundus) - R2 =", GRC_r2_2dp_NOPROF)),
                  (paste("Brienz/Thun (exl. profundus) - R2 =", GRC_B_r2_2dp)), 
                  (paste("Lucerne - R2 =", GRC_L_r2_2dp)), 
                  (paste("Constance - R2 =", GRC_C_r2_2dp)), 
                  (paste("Neuchatel/Biel - R2 =", GRC_N_r2_2dp)), 
                  (paste("Walen/Zurich - R2 =", GRC_B_r2_2dp))),
       col = c("black",
               "grey",
               "chartreuse4", 
               "#ff4489", 
               "goldenrod", 
               "turquoise4", 
               "purple4"),
       pt.cex = 1.2, cex = 0.8, pch = "-")
#dev.off()

#standard length similar 
par(mar=c(5,4,4,1)+.1)

SL <- pca

SL$SL_plot <- c()
for (i in 1:nrow(SL)){
  indiv_name <- as.character(SL$ind)[i]
  SL_sub <- subset(background, as.character(background$Indiv) == indiv_name)
  SL$SL_plot[i] <- as.character(SL_sub$standard_length)
}

SL <- subset(SL, SL$SL_plot != "missing")
SL$SL_plot <- as.numeric(SL$SL_plot)

PC1_SL <- lm(SL$SL_plot ~ SL$PC1)
summary(PC1_SL)
#Multiple R-squared:  0.3846,	Adjusted R-squared:  0.3777 
#F-statistic: 55.63 on 1 and 89 DF,  p-value: 5.527e-11

summary(PC1_SL)$r.squared

SL_r2 <- summary(PC1_SL)$r.squared
SL_r2_2dp <- format(round(SL_r2, 3), nsmall = 3)

plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "standard length",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. standard length - R2 =", SL_r2_2dp))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 3)

####Luzern
luzern_SL <- subset(SL, SL$loc == "Lucerne")
par(mar=c(5,4,4,1)+.1)

PC1_SL_L <- lm(luzern_SL$SL_plot ~ luzern_SL$PC1)
summary(PC1_SL_L)
#Multiple R-squared:  0.4546,	Adjusted R-squared:  0.4206 
#F-statistic: 13.34 on 1 and 16 DF,  p-value: 0.002149

SL_L_r2 <- summary(PC1_SL_L)$r.squared
SL_L_r2_2dp <- format(round(SL_L_r2, 3), nsmall = 3)


plot(luzern_SL$PC1, luzern_SL$SL_plot, col = as.character(luzern_SL$colour_plot), 
     pch=luzern_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("Luzern PC1 vs. SL - R2 =", SL_L_r2_2dp))
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")

####Walen/Zurich
walen_SL <- subset(SL, SL$loc == "Walen" | SL$loc == "Zurich")
PC1_SL_W <- lm(walen_SL$SL_plot ~ walen_SL$PC1)
summary(PC1_SL_W)
#Multiple R-squared:  0.6986,	Adjusted R-squared:  0.6798 
#F-statistic: 37.08 on 1 and 16 DF,  p-value: 1.564e-05

SL_W_r2 <- summary(PC1_SL_W)$r.squared
SL_W_r2_2dp <- format(round(SL_W_r2, 3), nsmall = 3)


plot(walen_SL$PC1, walen_SL$SL_plot, col = as.character(walen_SL$colour_plot), 
     pch=walen_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("walen PC1 vs. SL - R2 =", SL_W_r2_2dp))
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")


####Brienz/Thun - without profundal
brienz_SL <- subset(SL, SL$loc == "Brienz" | SL$loc == "Thun")
PC1_SL_B <- lm(brienz_SL$SL_plot ~ brienz_SL$PC1)
summary(PC1_SL_B)
#Multiple R-squared:  0.3913,	Adjusted R-squared:  0.3717 
#F-statistic: 19.93 on 1 and 31 DF,  p-value: 9.903e-05

SL_B_r2 <- summary(PC1_SL_B)$r.squared
SL_B_r2_2dp <- format(round(SL_B_r2, 3), nsmall = 3)


plot(brienz_SL$PC1, brienz_SL$SL_plot, col = as.character(brienz_SL$colour_plot), 
     pch=brienz_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("brienz/thun PC1 vs. SL - R2 =", SL_B_r2_2dp))
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")

####Neuchatel/biel
neuchatel_SL <- subset(SL, SL$loc == "Neuenburg" | SL$loc == "Biel")
PC1_SL_N <- lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1)
summary(PC1_SL_N)
#Multiple R-squared:  0.6771,	Adjusted R-squared:  0.6448 
#F-statistic: 20.97 on 1 and 10 DF,  p-value: 0.001012

SL_N_r2 <- summary(PC1_SL_N)$r.squared
SL_N_r2_2dp <- format(round(SL_N_r2, 3), nsmall = 3)


plot(neuchatel_SL$PC1, neuchatel_SL$SL_plot, col = as.character(neuchatel_SL$colour_plot), 
     pch=neuchatel_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("neuchatel/thun PC1 vs. SL - R2 =", SL_N_r2_2dp))
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")


####constance
constance_SL <- subset(SL, SL$loc == "Constance")
PC1_SL_C <- lm(constance_SL$SL_plot ~ constance_SL$PC1)
summary(PC1_SL_C)
#Multiple R-squared:  0.0618,	Adjusted R-squared:  -0.05547 
#F-statistic: 0.527 on 1 and 8 DF,  p-value: 0.4886

SL_C_r2 <- summary(PC1_SL_C)$r.squared
SL_C_r2_2dp <- format(round(SL_C_r2, 3), nsmall = 3)


plot(constance_SL$PC1, constance_SL$SL_plot, col = as.character(constance_SL$colour_plot), 
     pch=constance_SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("constance/thun PC1 vs. SL - R2 =", SL_C_r2_2dp))
abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")


#now plot full plot again with all ab lines separately
#tiff("Parallel/Parallel_css/All_342PC1_vs_StandardLength_noscaff22.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. SL"))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(luzern_SL$SL_plot ~ luzern_SL$PC1), lty = 3, col = "#ff4489")
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")
abline(lm(brienz_SL$SL_plot ~ brienz_SL$PC1), lty = 3, col = "chartreuse4")
abline(lm(neuchatel_SL$SL_plot ~ neuchatel_SL$PC1), lty = 3, col = "turquoise4")
#abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")
legend('topleft', 
       legend = c((paste("Overall - R2 =", SL_r2_2dp)),
                  (paste("Brienz/Thun - R2 =", SL_B_r2_2dp)), 
                  (paste("Lucerne - R2 =", SL_L_r2_2dp)), 
                  (paste("Constance (non-significant) - R2 =", SL_C_r2_2dp)), 
                  (paste("Neuchatel/Biel - R2 =", SL_N_r2_2dp)), 
                  (paste("Walen/Zurich - R2 =", SL_B_r2_2dp))),
       col = c("black",
               "chartreuse4", 
               "#ff4489", 
               "goldenrod", 
               "turquoise4", 
               "purple4"),
       pt.cex = 1.2, cex = 0.8, pch = "-")
#dev.off()
