###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	6.  parallel
#   - 6.1 parallel PCA
#		- 6.2 Fst
#		- 6.3 CSS
### see other script for		- 6.4 outlier analysis

#set directory and load background files
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)
chroms <- read.csv(file = "Parallel/chrom_names.txt", header = FALSE)
chromlist <- as.list(as.character(chroms$V1))

################# 6. 1 parallel PCA ############################### 
#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca <- read_table2("Parallel/parallel_filt.eigenvec", col_names = FALSE)
eigenval <- scan("Parallel/parallel_filt.eigenval")

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
#tiff("PCAs/plots/99PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("PCAs/plots/99PC1PC2_names.tiff", height=6, width=6, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P2")
#text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
#dev.off()

#tiff("PCAs/plots/99PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("PCAs/plots/99PC1PC3_names.tiff", height=6, width=6, units="in", res=300, compression="lzw")
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
#dev.off()

#tiff("Parallel/ParallelPC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
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
#legend('topright', 
#       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
#       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
#       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
legend('topright', 
       legend = c("Brienz", "Walen", "Luzern", "Neuchatel", "albeli",  "balchen"), 
       col = c("chartreuse2", "mediumpurple1",  "#ff4489", "turquoise4", "black", "black"), 
       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 21, 22))

#dev.off()

luzern_HGR <- subset(pca, pca$spp_loc == "C.zugensis_Lucerne")
luzern_LGR <- subset(pca, pca$spp_loc == "C.bodenbalchen_Lucerne")
segments(mean(luzern_LGR$PC2), mean(luzern_LGR$PC3), mean(luzern_HGR$PC2), mean(luzern_HGR$PC3),
         col = "black", lty = 3, lwd = 1)

brienz_HGR <- subset(pca, pca$spp_loc == "C.albellus_Brienz")
brienz_LGR <- subset(pca, pca$spp_loc == "C.alpinus_Brienz")
segments(mean(brienz_LGR$PC2), mean(brienz_LGR$PC3), mean(brienz_HGR$PC2), mean(brienz_HGR$PC3),
         col = "black", lty = 3, lwd = 1)

neuchatel_HGR <- subset(pca, pca$spp_loc == "C.candidus_Neuenburg")
neuchatel_LGR <- subset(pca, pca$spp_loc == "C.palaea_Neuenburg")
segments(mean(neuchatel_LGR$PC2), mean(neuchatel_LGR$PC3), mean(neuchatel_HGR$PC2), mean(neuchatel_HGR$PC3),
         col = "black", lty = 3, lwd = 1)

walen_HGR <- subset(pca, pca$spp_loc == "C.heglingus_Walen")
walen_LGR <- subset(pca, pca$spp_loc == "C.duplex_Walen")
segments(mean(walen_LGR$PC2), mean(walen_LGR$PC3), mean(walen_HGR$PC2), mean(walen_HGR$PC3),
         col = "black", lty = 3, lwd = 1)
#dev.off()


# CALULTE SLOPE ANGLE FOR EACH HGR/LGR COMPARISON
#get gradient of the slope following: slope = (y2 - y1)/(x2 - x1)
slope.luzern <- (mean(luzern_LGR$PC3)-mean(luzern_HGR$PC3))/(mean(luzern_LGR$PC2) - mean(luzern_HGR$PC2))
#then to get angle with x:
#angle theta with the x-axis tan(theta) = slope = (change in y) / (change in x)
#theta = tan_inverse (slope)
theta1 = atan(slope.luzern)

slope.brienz <- (mean(brienz_LGR$PC3)-mean(brienz_HGR$PC3))/(mean(brienz_LGR$PC2) - mean(brienz_HGR$PC2))
theta2 = atan(slope.brienz)

slope.walen <- (mean(walen_LGR$PC3)-mean(walen_HGR$PC3))/(mean(walen_LGR$PC2) - mean(walen_HGR$PC2))
theta3 = atan(slope.walen)

slope.neuchatel <- (mean(neuchatel_LGR$PC3)-mean(neuchatel_HGR$PC3))/(mean(neuchatel_LGR$PC2) - mean(neuchatel_HGR$PC2))
theta4 = atan(slope.neuchatel)


################# 6. 2 parallel fst ############################### 

#parallel fst outliers
brienz <- read.csv("Parallel/C.albellus_Brienz.indivlist_C.alpinus_Brienz.indivlist_fst.windowed.weir.fst", sep = "\t")
luzern <- read.csv("Parallel/C.zugensis_Lucerne.indivlist_C.bodenbalchen_Lucerne.indivlist_fst.windowed.weir.fst", sep = "\t")
neuchatel <- read.csv("Parallel/C.candidus_Neuenburg.indivlist_C.palaea_Neuenburg.indivlist_fst.windowed.weir.fst", sep = "\t")
walen <- read.csv("Parallel/C.heglingus_Walen.indivlist_C.duplex_Walen.indivlist_fst.windowed.weir.fst", sep = "\t")

brienz <- subset(brienz, as.integer(brienz$N_VARIANTS) > 10)
luzern <- subset(luzern, as.integer(luzern$N_VARIANTS) > 10)
neuchatel <- subset(neuchatel, as.integer(neuchatel$N_VARIANTS) > 10)
walen <- subset(walen, as.integer(walen$N_VARIANTS) > 10)

#par(mfrow=c(10,4))
#par(mar=c(1,1,1,1))
windowsize <- 50000
palette(c("#008958", "#ff4489", "#81b300", "#713c89"))
for (i in 1:40){
  chromosome <- as.character(chromlist[i])
  check_chrom <- levels(brienz$CHROM)
  check_chrom <- as.character(check_chrom)
  if (as.character(chromosome) %in% check_chrom){
    brienz_chrom <- subset(brienz, as.character(brienz$CHROM) == as.character(chromosome))
    brienz_chrom$lake <- "brienz"
    luzern_chrom <- subset(luzern, as.character(luzern$CHROM) == as.character(chromosome))  
    luzern_chrom$lake <- "luzern"
    neuchatel_chrom <- subset(neuchatel, as.character(neuchatel$CHROM) == as.character(chromosome))
    neuchatel_chrom$lake <- "neuchatel"
    walen_chrom <- subset(walen, as.character(walen$CHROM) == as.character(chromosome))
    walen_chrom$lake <- "walen"
    fulldf <- rbind(brienz_chrom, luzern_chrom, neuchatel_chrom, walen_chrom)
    fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
    fulldf$lake <- as.factor(fulldf$lake)
    fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
    #shapes = c(19, 17, 18, 16) 
    #shapes <- shapes[as.numeric(fulldf$lake)]
    #plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$lake), main = chromosome, ylim = c(0,1), pch = shapes)
    options(scipen=8)
    plotname <- paste("CHR", as.character(i), "_", as.character(chromosome), sep = "")
    #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
    plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$lake), main = plotname, ylim = c(0,1))
    legend("topright", legend = levels(fulldf$lake),
           #col =  c("black", "red", "green", "blue"),
           col = c("#008958", "#ff4489", "#81b300", "#713c89"),
           pch = c(16))
    #dev.off()
  }
}

palette(c("black", "red"))
sharedoutliers_total = data.frame()
for (i in 1:40){
  chromosome <- as.character(chromlist[i])
  check_chrom <- levels(brienz$CHROM)
  check_chrom <- as.character(check_chrom)
  if (as.character(chromosome) %in% check_chrom){
    brienz_chrom <- subset(brienz, as.character(brienz$CHROM) == as.character(chromosome))
    brienz_chrom$lake <- "brienz"
    luzern_chrom <- subset(luzern, as.character(luzern$CHROM) == as.character(chromosome))  
    luzern_chrom$lake <- "luzern"
    neuchatel_chrom <- subset(neuchatel, as.character(neuchatel$CHROM) == as.character(chromosome))
    neuchatel_chrom$lake <- "neuchatel"
    walen_chrom <- subset(walen, as.character(walen$CHROM) == as.character(chromosome))
    walen_chrom$lake <- "walen"
    fulldf <- rbind(brienz_chrom, luzern_chrom, neuchatel_chrom, walen_chrom)
    fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
    fulldf$lake <- as.factor(fulldf$lake)
    fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
    fulldf$shared <- "N"
    for (j in 1:nrow(fulldf)){
      window_interest <- (fulldf$averagebp)[j]
      shared_subset <- subset(fulldf, as.character(fulldf$averagebp) == as.character(window_interest))
      shared_subset_fstcutoff <- subset(shared_subset, shared_subset$MEAN_FST > 0.6)
      if (nrow(shared_subset_fstcutoff) > 1){
        fulldf$shared[j] <- "Y"
      }
    }
    fulldf$shared <- as.factor(fulldf$shared)
    shared_outliers <- subset(fulldf, as.character(fulldf$shared) == "Y")
    sharedoutliers_total <- rbind(sharedoutliers_total, shared_outliers)
    shapes = c(1, 16) 
    shapes <- shapes[as.numeric(fulldf$shared)]
    options(scipen=8)
    plotname <- paste("CHR", as.character(i), "_", as.character(chromosome), "_highlighted0.6outliers", sep = "")
    #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
    plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$shared), pch = shapes, main = plotname, ylim = c(0,1), xlab = "bp", ylab = "mean Fst")
    legend("topright", legend = levels(fulldf$shared),
           col =  c("black", "red"),
           pch = c(1, 16))
    #dev.off()
  }
}

sharedoutliers_total_brienz <- subset(sharedoutliers_total, as.character(sharedoutliers_total$lake) == "brienz")
#write.csv(sharedoutliers_total_brienz,"sharedoutliers.csv", row.names = FALSE)

# Saving on object in RData format
save(sharedoutliers_total, file = "sharedoutliers_total.RData")

# To load the data again
#load("sharedoutliers_total.RData")


#and the same including pooled samples
####################################################################################
#
pooled <- read.csv("Parallel/HGR.indivlist_LGR.indivlist_fst.windowed.weir.fst", sep = "\t")
brienz <- read.csv("Parallel/C.albellus_Brienz.indivlist_C.alpinus_Brienz.indivlist_fst.windowed.weir.fst", sep = "\t")
luzern <- read.csv("Parallel/C.zugensis_Lucerne.indivlist_C.bodenbalchen_Lucerne.indivlist_fst.windowed.weir.fst", sep = "\t")
neuchatel <- read.csv("Parallel/C.candidus_Neuenburg.indivlist_C.palaea_Neuenburg.indivlist_fst.windowed.weir.fst", sep = "\t")
walen <- read.csv("Parallel/C.heglingus_Walen.indivlist_C.duplex_Walen.indivlist_fst.windowed.weir.fst", sep = "\t")

pooled <- subset(pooled, as.integer(pooled$N_VARIANTS) > 10)
brienz <- subset(brienz, as.integer(brienz$N_VARIANTS) > 10)
luzern <- subset(luzern, as.integer(luzern$N_VARIANTS) > 10)
neuchatel <- subset(neuchatel, as.integer(neuchatel$N_VARIANTS) > 10)
walen <- subset(walen, as.integer(walen$N_VARIANTS) > 10)

#par(mfrow=c(10,4))
#par(mar=c(1,1,1,1))
windowsize <- 50000
palette(c("#008958", "#ff4489", "#81b300", "black", "#713c89"))
for (i in 1:40){
  chromosome <- as.character(chromlist[i])
  check_chrom <- levels(brienz$CHROM)
  check_chrom <- as.character(check_chrom)
  if (as.character(chromosome) %in% check_chrom){
    brienz_chrom <- subset(brienz, as.character(brienz$CHROM) == as.character(chromosome))
    brienz_chrom$lake <- "brienz"
    luzern_chrom <- subset(luzern, as.character(luzern$CHROM) == as.character(chromosome))  
    luzern_chrom$lake <- "luzern"
    neuchatel_chrom <- subset(neuchatel, as.character(neuchatel$CHROM) == as.character(chromosome))
    neuchatel_chrom$lake <- "neuchatel"
    walen_chrom <- subset(walen, as.character(walen$CHROM) == as.character(chromosome))
    walen_chrom$lake <- "walen"
    pooled_chrom <- subset(pooled, as.character(pooled$CHROM) == as.character(chromosome))
    pooled_chrom$lake <- "pooled"
    fulldf <- rbind(brienz_chrom, luzern_chrom, neuchatel_chrom, walen_chrom, pooled_chrom)
    fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
    fulldf$lake <- as.factor(fulldf$lake)
    fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
    #shapes = c(19, 17, 18, 16) 
    #shapes <- shapes[as.numeric(fulldf$lake)]
    #plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$lake), main = chromosome, ylim = c(0,1), pch = shapes)
    plotname <- paste("CHR", as.character(i), "_", as.character(chromosome), sep = "")
    #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
    plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$lake), main = plotname, ylim = c(0,1))
    legend("topright", legend = levels(fulldf$lake),
           #col =  c("black", "red", "green", "blue"),
           col = c("#008958", "#ff4489", "#81b300", "black", "#713c89"),
           pch = c(16))
    #dev.off()
  }
}

palette(c("black", "red"))
sharedoutliers_total_pooled = data.frame()
for (i in 1:40){
  chromosome <- as.character(chromlist[i])
  check_chrom <- levels(brienz$CHROM)
  check_chrom <- as.character(check_chrom)
  if (as.character(chromosome) %in% check_chrom){
    brienz_chrom <- subset(brienz, as.character(brienz$CHROM) == as.character(chromosome))
    brienz_chrom$lake <- "brienz"
    luzern_chrom <- subset(luzern, as.character(luzern$CHROM) == as.character(chromosome))  
    luzern_chrom$lake <- "luzern"
    neuchatel_chrom <- subset(neuchatel, as.character(neuchatel$CHROM) == as.character(chromosome))
    neuchatel_chrom$lake <- "neuchatel"
    walen_chrom <- subset(walen, as.character(walen$CHROM) == as.character(chromosome))
    walen_chrom$lake <- "walen"
    pooled_chrom <- subset(pooled, as.character(pooled$CHROM) == as.character(chromosome))
    pooled_chrom$lake <- "pooled"
    fulldf <- rbind(brienz_chrom, luzern_chrom, neuchatel_chrom, walen_chrom, pooled_chrom)
    fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
    fulldf$lake <- as.factor(fulldf$lake)
    fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
    fulldf$shared <- "N"
    for (j in 1:nrow(fulldf)){
      window_interest <- (fulldf$averagebp)[j]
      shared_subset <- subset(fulldf, as.character(fulldf$averagebp) == as.character(window_interest))
      shared_subset_fstcutoff <- subset(shared_subset, shared_subset$MEAN_FST > 0.6)
      if (nrow(shared_subset_fstcutoff) > 2){
        fulldf$shared[j] <- "Y"
      }
    }
    fulldf$shared <- as.factor(fulldf$shared)
    shared_outliers <- subset(fulldf, as.character(fulldf$shared) == "Y")
    sharedoutliers_total_pooled <- rbind(sharedoutliers_total_pooled, shared_outliers)
    shapes = c(1, 16) 
    shapes <- shapes[as.numeric(fulldf$shared)]
    options(scipen=8)
    plotname <- paste("CHR", as.character(i), "_", as.character(chromosome), "_highlighted0.6outliers", sep = "")
    #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
    plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$shared), pch = shapes, main = plotname, ylim = c(0,1), xlab = "bp", ylab = "mean Fst")
    legend("topright", legend = levels(fulldf$shared),
           col =  c("black", "red"),
           pch = c(1, 16))
    #dev.off()
  }
}

#now produce a file where the pooled fst resolves a pattern observed in the other sampless
sharedoutliers_total_ONLY_POOLEDandSHARED <- subset(sharedoutliers_total_pooled, as.character(sharedoutliers_total_pooled$lake) == "pooled")
#write.csv(sharedoutliers_total_ONLY_POOLEDandSHARED,"sharedoutliers_pooled.csv", row.names = FALSE)

#just plotting the pooled sample Fsts
for (i in 1:40){
  chromosome <- as.character(chromlist[i])
  check_chrom <- levels(brienz$CHROM)
  check_chrom <- as.character(check_chrom)
  if (as.character(chromosome) %in% check_chrom){
    pooled_chrom <- subset(pooled, as.character(pooled$CHROM) == as.character(chromosome))
    pooled_chrom$lake <- "pooled"
    fulldf <- rbind(pooled_chrom)
    fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
    fulldf$lake <- as.factor(fulldf$lake)
    fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
    #shapes = c(19, 17, 18, 16) 
    #shapes <- shapes[as.numeric(fulldf$lake)]
    #plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$lake), main = chromosome, ylim = c(0,1), pch = shapes)
    plotname <- paste("CHR", as.character(i), "pooled", sep = "")
    #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
    plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = "black", main = plotname, ylim = c(0,1))
    abline(h = 0.6, col = "red")
    #dev.off()
  }
}


########################################################
#   STATS OF OUTLIERS
#get mean snps in windows for full data and outliers
par(mfrow=c(2,2))
brienz$fst_mean <- brienz$MEAN_FST
brienz$fst_mean[brienz$fst_mean < 0] <- 0
hist(brienz$fst_mean, xlim = c(0,1), ylim = c(0,25000), xlab = "Fst", main = "brienz fst")

walen$fst_mean <- walen$MEAN_FST
walen$fst_mean[walen$fst_mean < 0] <- 0
hist(walen$fst_mean, xlim = c(0,1), ylim = c(0,25000), xlab = "Fst", main = "walen fst")

neuchatel$fst_mean <- neuchatel$MEAN_FST
neuchatel$fst_mean[neuchatel$fst_mean < 0] <- 0
hist(neuchatel$fst_mean, xlim = c(0,1), ylim = c(0,25000), xlab = "Fst", main = "neuchatel fst")


luzern$fst_mean <- luzern$MEAN_FST
luzern$fst_mean[luzern$fst_mean < 0] <- 0
hist(luzern$fst_mean, xlim = c(0,1), ylim = c(0,25000), xlab = "Fst", main = "luzern fst")


##BRIENZ
mean_snps_brienz <- mean(brienz$N_VARIANTS)
mean_fst_brienz <- mean(brienz$fst_mean)
#mean snps full
mean_snps_brienz
mean_fst_brienz
brienz_outliers <- subset(sharedoutliers_total, as.character(sharedoutliers_total$lake) == "brienz")
mean_snps_brienz_outliers <- mean(brienz_outliers$N_VARIANTS)
#mean in outliers
mean_snps_brienz_outliers
#max
max(brienz_outliers$N_VARIANTS)
#min
min(brienz_outliers$N_VARIANTS)
#mean outlier fst
mean(brienz_outliers$MEAN_FST)


##WALEN
mean_snps_walen <- mean(walen$N_VARIANTS)
mean_fst_walen <- mean(walen$fst_mean)
#mean snps full
mean_snps_walen
mean_fst_walen
walen_outliers <- subset(sharedoutliers_total, as.character(sharedoutliers_total$lake) == "walen")
mean_snps_walen_outliers <- mean(walen_outliers$N_VARIANTS)
#mean in outliers
mean_snps_walen_outliers
#max
max(walen_outliers$N_VARIANTS)
#min
min(walen_outliers$N_VARIANTS)
#mean outlier fst
mean(walen_outliers$MEAN_FST)


##NEUCHATEL
mean_snps_neuchatel <- mean(neuchatel$N_VARIANTS)
mean_fst_neuchatel <- mean(neuchatel$fst_mean)
#mean snps full
mean_snps_neuchatel
mean_fst_neuchatel
neuchatel_outliers <- subset(sharedoutliers_total, as.character(sharedoutliers_total$lake) == "neuchatel")
mean_snps_neuchatel_outliers <- mean(neuchatel_outliers$N_VARIANTS)
#mean in outliers
mean_snps_neuchatel_outliers
#max
max(neuchatel_outliers$N_VARIANTS)
#min
min(neuchatel_outliers$N_VARIANTS)
#mean outlier fst
mean(neuchatel_outliers$MEAN_FST)


##LUZERN
mean_snps_luzern <- mean(luzern$N_VARIANTS)
mean_fst_luzern <- mean(luzern$fst_mean)
#mean snps full
mean_snps_luzern
mean_fst_luzern
luzern_outliers <- subset(sharedoutliers_total, as.character(sharedoutliers_total$lake) == "luzern")
mean_snps_luzern_outliers <- mean(luzern_outliers$N_VARIANTS)
#mean in outliers
mean_snps_luzern_outliers
#max
max(luzern_outliers$N_VARIANTS)
#min
min(luzern_outliers$N_VARIANTS)
#mean outlier fst
mean(luzern_outliers$MEAN_FST)

par(mfrow=c(2,2))
hist(brienz_outliers$MEAN_FST, xlim = c(0,1), ylim = c(0,25), xlab = "Fst", main = "brienz shared outliers fst")
hist(walen_outliers$MEAN_FST, xlim = c(0,1), ylim = c(0,25), xlab = "Fst", main = "walen shared outliers fst")
hist(neuchatel_outliers$MEAN_FST, xlim = c(0,1), ylim = c(0,25), xlab = "Fst", main = "neuchatel shared outliers fst")
hist(luzern_outliers$MEAN_FST, xlim = c(0,1), ylim = c(0,25), xlab = "Fst", main = "luzern shared outliers fst")

################################################################
#finding fixed differences:
brienz_fixed <- subset(brienz, brienz$MEAN_FST == 1)
print("brienz")
print(nrow(brienz_fixed))
walen_fixed <- subset(walen, walen$MEAN_FST == 1)
print("walen")
print(nrow(walen_fixed))
neuchatel_fixed <- subset(neuchatel, neuchatel$MEAN_FST == 1)
print("neuchatel")
print(nrow(neuchatel_fixed))
luzern_fixed <- subset(luzern, luzern$MEAN_FST == 1)
print("luzern")
print(nrow(luzern_fixed))


#       CALCULATE OUTLIERS - top 1% FOR EACH LAKE
#we want to calc outliers for each plot:
#brienz
#calc empirical cumulative distribution function of Fst within each lake
brienz_9 <- brienz[brienz$N_VARIANTS>9,]
FnE_brienz<-ecdf(brienz_9$MEAN_FST)
#then make new df with old data and the new p-value for each fst
brienz_OUT<-cbind(brienz_9,FST_p=(1-FnE_brienz(brienz_9$MEAN_FST)))
#then filter to keep only those with p<0.01
OUT_brienz<-brienz_OUT[brienz_OUT$FST_p<0.01,]
OUT_brienz

#luzern
luzern_9 <- luzern[luzern$N_VARIANTS>9,]
FnE_luzern<-ecdf(luzern_9$MEAN_FST)
luzern_OUT<-cbind(luzern_9,FST_p=(1-FnE_luzern(luzern_9$MEAN_FST)))
OUT_luzern<-luzern_OUT[luzern_OUT$FST_p<0.01,]
OUT_luzern

#walen
walen_9 <- walen[walen$N_VARIANTS>9,]
FnE_walen<-ecdf(walen_9$MEAN_FST)
walen_OUT<-cbind(walen_9,FST_p=(1-FnE_walen(walen_9$MEAN_FST)))
OUT_walen<-walen_OUT[walen_OUT$FST_p<0.01,]
OUT_walen

#neuchatel
neuchatel_9 <- neuchatel[neuchatel$N_VARIANTS>9,]
FnE_neuchatel<-ecdf(neuchatel_9$MEAN_FST)
neuchatel_OUT<-cbind(neuchatel_9,FST_p=(1-FnE_neuchatel(neuchatel_9$MEAN_FST)))
OUT_neuchatel<-neuchatel_OUT[neuchatel_OUT$FST_p<0.01,]
OUT_neuchatel

#now save outliers from each comparison:
write.csv(OUT_brienz, "Parallel/Parallel_fst/Brienz_outliers.csv")
write.csv(OUT_luzern, "Parallel/Parallel_fst/Lucerne_outliers.csv")
write.csv(OUT_walen, "Parallel/Parallel_fst/Walen_outliers.csv")
write.csv(OUT_neuchatel, "Parallel/Parallel_fst/Neuchatel_outliers.csv")


#write loop to go chromosome by chromosome, subset the outliers for each lake by chromomsome, rbind together then go window by window and say whether window is shared between 1, 2, 3, 4 lakes by woking out the row number
#or just cbind brien_OUT etc. and then add column for each lake denoting where outlier window is i.e. p< value and then count up the yesses in those columns

# outliers_total = data.frame()
# for (i in 1:40){
#   chromosome <- as.character(chromlist[i])
#   check_chrom <- levels(brienz$CHROM)
#   check_chrom <- as.character(check_chrom)
#   if (as.character(chromosome) %in% check_chrom){
#     brienz_chrom <- subset(OUT_brienz, as.character(OUT_brienz$CHROM) == as.character(chromosome))
#     brienz_chrom$lake <- "brienz"
#     luzern_chrom <- subset(OUT_luzern, as.character(OUT_luzern$CHROM) == as.character(chromosome))  
#     luzern_chrom$lake <- "luzern"
#     neuchatel_chrom <- subset(OUT_neuchatel, as.character(OUT_neuchatel$CHROM) == as.character(chromosome))
#     neuchatel_chrom$lake <- "neuchatel"
#     walen_chrom <- subset(OUT_walen, as.character(OUT_walen$CHROM) == as.character(chromosome))
#     walen_chrom$lake <- "walen"
#     fulldf <- rbind(brienz_chrom, luzern_chrom, neuchatel_chrom, walen_chrom)
#     fulldf$averagebp <- (fulldf$BIN_END-(windowsize/2))
#     fulldf$lake <- as.factor(fulldf$lake)
#     fulldf$MEAN_FST[fulldf$MEAN_FST < 0] <- 0
#     fulldf$shared <- "N"
#     for (j in 1:nrow(fulldf)){
#       window_interest <- (fulldf$averagebp)[j]
#       shared_subset <- subset(fulldf, as.character(fulldf$averagebp) == as.character(window_interest))
#       fulldf$shared[j] <- nrow(shared_subset)
#       }
#     fulldf$shared <- as.factor(fulldf$shared)
#     shared_outliers <- fulldf
#     outliers_total <- rbind(outliers_total, shared_outliers)
#     shapes = c(1, 16) 
#     shapes <- shapes[as.numeric(fulldf$shared)]
#     options(scipen=8)
#     plotname <- paste("CHR", as.character(i), "_", as.character(chromosome), "_highlighted0.6outliers", sep = "")
#     #tiff(paste(as.character(plotname), ".tiff", sep = ""), height=8, width=20, units="in", res=300, compression="lzw")
#     plot(x = fulldf$averagebp, y = fulldf$MEAN_FST, col = as.factor(fulldf$shared), pch = shapes, main = plotname, ylim = c(0,1), xlab = "bp", ylab = "mean Fst")
#     legend("topright", legend = levels(fulldf$shared),
#            col =  c("black", "red"),
#            pch = c(1, 16))
#     #dev.off()
#   }
# }

#load css data
css<- read.csv("Parallel/all_CSS_output_500kpermutations_50000basepair50000step.window.pca.full_noheaders.txt", sep = "\t", head = FALSE)

css_colnames6 <- colnames(pooled)
css_colnames5 <- css_colnames6[-c(5)]
colnames(css) <- c(css_colnames5, "q_value")

#now plot css scores vs q value
plot(css$MEAN_FST, css$q_value)

#css
css_23 <- css[css$N_VARIANTS>23,]
FnE_css<-ecdf(css_23$MEAN_FST)
css_OUT<-cbind(css_23,FST_p=(1-FnE_css(css_23$MEAN_FST)))
OUT_css<-css_OUT[css_OUT$FST_p<0.01,]
OUT_css

#and now filter top one percent by css score
#CSS_permuted_outliers <- subset(css, (1-(css$q_value)) < 0.005)
CSS_permuted_outliers <- subset(OUT_css, (1-(OUT_css$q_value)) < 0.01)

#write csv of css outliers
#OUTPUT
write.csv(CSS_permuted_outliers,"Parallel/CSS_1percent_outliers.csv", row.names = FALSE)

#find adjacent rows where >= 3 adjacent windows are outliers
test_adj <- CSS_permuted_outliers
test_adj$adj_3 <- "NO"
for (i in 1:((length(test_adj$adj_3))-2)){
  j <- as.numeric(i+1)
  k <- as.numeric(j+1)
  endi <- as.integer((test_adj$BIN_END[i]+1))
  startj <- as.integer(test_adj$BIN_START[j])
  endj <- as.integer((test_adj$BIN_END[j]+1))
  startk <- as.integer(test_adj$BIN_START[k])
  if (endi == startj){
    if (endj ==startk){
      test_adj$adj_3[i] <- "YES"
      test_adj$adj_3[j] <- "YES"
      test_adj$adj_3[k] <- "YES"}
  }
}

test_adj_filt <- subset(test_adj, test_adj$adj_3 == "YES")


#############
#     PLOT LANDSCAPE OUTLEIER PLOT SMALL POINTS
#now produce chromosome list to we can plot each chromosome
b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

#pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj.pdf",width=22, height=12)
pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj_fixed.pdf",width=22, height=12)
par(mfrow=c(5,1),mar=c(0.5,0.15,0,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
l <- as.vector(lengthsdf$V1)

#fix lengths for chromosomes that have partial duplicates where regions have been removed:
#PGA_scaffold3__454_contigs__length_92224161     1       31500000
#PGA_scaffold6__535_contigs__length_65391737     39000000        65391737
#PGA_scaffold16__334_contigs__length_54216998    41500000        54216998
#PGA_scaffold36__136_contigs__length_43663377    32500000        43663377
lnew <- l
lnew[4] <- l[4]-31500000

lnew[7] <- l[7]-(65391737-39000000)

lnew[17] <- l[17]-(54216998-41500000)

lnew[37] <- l[37]-(43663377-32500000)

#remove the absent chromosomes
lnew2 <- lnew[-c(22, 28, 32, 35, 38, 40)]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:170),nrow=5,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("slategray1", "slategray2")
altcols <- rep(altcols, as.integer(length(a)/2))

#now plot for each lake
#brienz
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- brienz_9$BIN_END[brienz_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- brienz_9$MEAN_FST[brienz_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="Fst - BR",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.5)
  #and add p-value outliers i.e. top 1% of fst values which we replot in pink
  outpoints_xFst <- OUT_brienz$BIN_END[OUT_brienz$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_brienz$MEAN_FST[OUT_brienz$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="darkgreen",pch=20,cex=1.5)
  if (i<2){
    axis(side=2)
  }
}

#luzern
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- luzern_9$BIN_END[luzern_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- luzern_9$MEAN_FST[luzern_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="Fst - LU",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_luzern$BIN_END[OUT_luzern$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_luzern$MEAN_FST[OUT_luzern$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="#ff4489",pch=20,cex=1.5)
  if (i<2){
    axis(side=2)
  }
}

#walen
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- walen_9$BIN_END[walen_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- walen_9$MEAN_FST[walen_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="Fst - WA",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_walen$BIN_END[OUT_walen$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_walen$MEAN_FST[OUT_walen$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="#713c89",pch=20,cex=1.5)
  if (i<2){
    axis(side=2)
  }
}

#neuchatel
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- neuchatel_9$BIN_END[neuchatel_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- neuchatel_9$MEAN_FST[neuchatel_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="Fst - NE",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_neuchatel$BIN_END[OUT_neuchatel$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_neuchatel$MEAN_FST[OUT_neuchatel$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="deepskyblue3",pch=20,cex=1.5)
  if (i<2){
    axis(side=2)
  }
}


altcols_pooled <- c("grey68", "grey50")
altcols_pooled <- rep(altcols_pooled, as.integer(length(a)/2))

#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css_23$BIN_END[css_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css_23$MEAN_FST[css_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.3),col=altcols_pooled[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_css$BIN_END[OUT_css$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css$MEAN_FST[OUT_css$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=1.5)
  outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
dev.off()

#############
#     PLOT LANDSCAPE OUTLEIER PLOT BIG POINTS
#now produce chromosome list to we can plot each chromosome
b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_big_fixed.pdf",width=10, height=6)

par(mfrow=c(5,1),mar=c(0.5,0.15,0,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
l <- as.vector(lengthsdf$V1)

#fix lengths for chromosomes that have partial duplicates where regions have been removed:
#PGA_scaffold3__454_contigs__length_92224161     1       31500000
#PGA_scaffold6__535_contigs__length_65391737     39000000        65391737
#PGA_scaffold16__334_contigs__length_54216998    41500000        54216998
#PGA_scaffold36__136_contigs__length_43663377    32500000        43663377
lnew <- l
lnew[4] <- l[4]-31500000

lnew[7] <- l[7]-(65391737-39000000)

lnew[17] <- l[17]-(54216998-41500000)

lnew[37] <- l[37]-(43663377-32500000)

#remove the absent chromosomes
lnew2 <- lnew[-c(22, 28, 32, 35, 38, 40)]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:170),nrow=5,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("slategray1", "slategray2")
altcols <- rep(altcols, as.integer(length(a)/2))

#now plot for each lake
#brienz
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- brienz_9$BIN_END[brienz_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- brienz_9$MEAN_FST[brienz_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=.5,ylab="Fst - BR",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.3)
  #and add p-value outliers i.e. top 1% of fst values which we replot in pink
  outpoints_xFst <- OUT_brienz$BIN_END[OUT_brienz$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_brienz$MEAN_FST[OUT_brienz$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="darkgreen",pch=20,cex=.5)
  if (i<2){
    axis(side=2, cex.axis = 0.5)
  }
}

#luzern
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- luzern_9$BIN_END[luzern_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- luzern_9$MEAN_FST[luzern_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=.5,ylab="Fst - LU",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.3)
  outpoints_xFst <- OUT_luzern$BIN_END[OUT_luzern$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_luzern$MEAN_FST[OUT_luzern$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="#ff4489",pch=20,cex=.5)
  if (i<2){
    axis(side=2, cex.axis = 0.5)
  }
}

#walen
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- walen_9$BIN_END[walen_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- walen_9$MEAN_FST[walen_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=.5,ylab="Fst - WA",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.3)
  outpoints_xFst <- OUT_walen$BIN_END[OUT_walen$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_walen$MEAN_FST[OUT_walen$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="#713c89",pch=20,cex=.5)
  if (i<2){
    axis(side=2, cex.axis = 0.5)
  }
}

#neuchatel
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- neuchatel_9$BIN_END[neuchatel_9$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- neuchatel_9$MEAN_FST[neuchatel_9$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=.5,ylab="Fst - NE",xaxt="n",axes=F,ylim=c(0,1),col=altcols[i], cex.axis = 0.3)
  outpoints_xFst <- OUT_neuchatel$BIN_END[OUT_neuchatel$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_neuchatel$MEAN_FST[OUT_neuchatel$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col="deepskyblue3",pch=20,cex=.5)
  if (i<2){
    axis(side=2, cex.axis = 0.5)
  }
}


altcols_pooled <- c("grey80", "grey60")
altcols_pooled <- rep(altcols_pooled, as.integer(length(a)/2))

#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css_23$BIN_END[css_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css_23$MEAN_FST[css_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=.5,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.3),col=altcols_pooled[i], cex.axis = 0.3)
  outpoints_xFst <- OUT_css$BIN_END[OUT_css$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css$MEAN_FST[OUT_css$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=.5)
  outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=.5)
  if (i<2){
    axis(side=2, cex.axis = 0.5)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col="snow4")
}
dev.off()

#to get r2 between fst of windows between lakes:
#walen has the most windows present so use this to full_join the others:
library(dplyr)
walen_pre <- walen
walen_pre$loc <- paste(walen_pre$CHROM, walen_pre$BIN_START)
walen_pre <-  walen_pre[,7:8]
colnames(walen_pre) <- c("walen_fst", "loc")
luzern_pre <- luzern
luzern_pre$loc <- paste(luzern_pre$CHROM, luzern_pre$BIN_START)
luzern_pre <-  luzern_pre[,7:8]
colnames(luzern_pre) <- c("luzern_fst", "loc")
brienz_pre <- brienz
brienz_pre$loc <- paste(brienz_pre$CHROM, brienz_pre$BIN_START)
brienz_pre <-  brienz_pre[,7:8]
colnames(brienz_pre) <- c("brienz_fst", "loc")
neuchatel_pre <- neuchatel
neuchatel_pre$loc <- paste(neuchatel_pre$CHROM, neuchatel_pre$BIN_START)
neuchatel_pre <- neuchatel_pre[,7:8]
colnames(neuchatel_pre) <- c("neuchatel_fst", "loc")


w_l <-  as.data.frame(walen_pre$loc)
w_l$walen_fst <- walen$fst_mean
colnames(w_l) <- c("loc", "walen_fst")
w_l$loc <- as.character(w_l$loc)

test <- full_join(w_l, luzern_pre, by = "loc")
test2<- full_join(test, brienz_pre, by = "loc")
test3 <- full_join(test2, neuchatel_pre, by = "loc")

test4 <- na.omit(test3)

plot(test4$walen_fst, test4$luzern_fst)
w_l.lm = lm(test4$walen_fst ~ test4$luzern_fst)
summary(w_l.lm)
summary(w_l.lm)$r.squared 

plot(test4$walen_fst, test4$brienz_fst)
w_b.lm = lm(test4$walen_fst ~ test4$brienz_fst)
summary(w_b.lm)
summary(w_b.lm)$r.squared

plot(test4$walen_fst, test4$neuchatel_fst)
w_n.lm = lm(test4$walen_fst ~ test4$neuchatel_fst)
summary(w_n.lm)
summary(w_n.lm)$r.squared

plot(test4$luzern_fst, test4$brienz_fst)
l_b.lm = lm(test4$luzern_fst ~ test4$brienz_fst)
summary(l_b.lm)
summary(l_b.lm)$r.squared

plot(test4$luzern_fst, test4$neuchatel_fst)
l_n.lm = lm(test4$luzern_fst ~ test4$neuchatel_fst)
summary(l_n.lm)
summary(l_n.lm)$r.squared

plot(test4$brienz_fst, test4$neuchatel_fst)
b_n.lm = lm(test4$brienz_fst ~ test4$neuchatel_fst)
summary(b_n.lm)
summary(b_n.lm)$r.squared

####and now the same but just for outliers
#to get r2 between fst of windows between lakes:
#walen has the most windows present so use this to full_join the others:
library(dplyr)
walen_pre <- OUT_walen
walen_pre$loc <- paste(walen_pre$CHROM, walen_pre$BIN_START)
walen_pre <-  walen_pre[,7:9]
walen_pre <-  walen_pre[,-2]
colnames(walen_pre) <- c("walen_fst", "loc")

luzern_pre <- OUT_luzern
luzern_pre$loc <- paste(luzern_pre$CHROM, luzern_pre$BIN_START)
luzern_pre <-  luzern_pre[,7:9]
luzern_pre <-  luzern_pre[,-2]
colnames(luzern_pre) <- c("luzern_fst", "loc")

brienz_pre <- OUT_brienz
brienz_pre$loc <- paste(brienz_pre$CHROM, brienz_pre$BIN_START)
brienz_pre <-  brienz_pre[,7:9]
brienz_pre <-  brienz_pre[,-2]
colnames(brienz_pre) <- c("brienz_fst", "loc")

neuchatel_pre <- OUT_neuchatel
neuchatel_pre$loc <- paste(neuchatel_pre$CHROM, neuchatel_pre$BIN_START)
neuchatel_pre <- neuchatel_pre[,7:9]
neuchatel_pre <- neuchatel_pre[,-2]
colnames(neuchatel_pre) <- c("neuchatel_fst", "loc")

w_l <-  as.data.frame(walen_pre$loc)
w_l$walen_fst <- OUT_walen$fst_mean
colnames(w_l) <- c("loc", "walen_fst")
w_l$loc <- as.character(w_l$loc)

test <- full_join(w_l, luzern_pre, by = "loc")
test2<- full_join(test, brienz_pre, by = "loc")
test3 <- full_join(test2, neuchatel_pre, by = "loc")

WL <- test3[,2:3]
WL$loc <- test3$loc
WL <- na.omit(WL)
plot(WL$walen_fst, WL$luzern_fst)
w_l.lm = lm(WL$luzern_fst~WL$walen_fst)
abline(lm(WL$luzern_fst~WL$walen_fst))
summary(w_l.lm)
summary(w_l.lm)$r.squared 

WB <- test3[,2:4]
WB <- WB[,-2]
WB$loc <- test3$loc
WB <- na.omit(WB)
plot(WB$walen_fst, WB$brienz_fst)
w_b.lm = lm(WB$brienz_fst~WB$walen_fst)
abline(lm(WB$brienz_fst~WB$walen_fst))
summary(w_b.lm)
summary(w_b.lm)$r.squared

WN <- test3[,2:5]
WN <- WN[,-2]
WN <- WN[,-2]
WN$loc <- test3$loc
WN <- na.omit(WN)
plot(WN$walen_fst, WN$neuchatel_fst)
w_n.lm = lm(WN$neuchatel_fst~WN$walen_fst)
abline(lm(WN$neuchatel_fst~WN$walen_fst))
summary(w_n.lm)
summary(w_n.lm)$r.squared

LB <- test3[,3:4]
LB$loc <- test3$loc
LB <- na.omit(LB)
plot(LB$luzern_fst, LB$brienz_fst)
l_b.lm = lm(LB$brienz_fst~LB$luzern_fst)
abline(lm(LB$brienz_fst~LB$luzern_fst))
summary(l_b.lm)
summary(l_b.lm)$r.squared

LN <- test3[,3:5]
LN <- LN[,-2]
LN$loc <- test3$loc
LN <- na.omit(LN)
plot(LN$luzern_fst, LN$neuchatel_fst)
l_n.lm = lm(LN$neuchatel_fst~LN$luzern_fst)
abline(lm(LN$neuchatel_fst~LN$luzern_fst))
summary(l_n.lm)
summary(l_n.lm)$r.squared

BN <- test3[,4:5]
BN$loc <- test3$loc
BN <- na.omit(BN)
plot(BN$brienz_fst, BN$neuchatel_fst)
b_n.lm = lm(BN$brienz_fst ~ BN$neuchatel_fst)
abline(lm(BN$neuchatel_fst~BN$brienz_fst))
summary(b_n.lm)
summary(b_n.lm)$r.squared

#get outlier windows which are shared between lakes
test3$missing <- rowSums(is.na(test3))
common2 <- subset(test3, test3$missing < 3)
common1 <- subset(test3, test3$missing < 2)

write.csv(common2, "shared_outliers_narrowsense.csv")

filled <- common2
walen$loc <- paste(walen$CHROM, walen$BIN_START)
luzern$loc <- paste(luzern$CHROM, luzern$BIN_START)
brienz$loc <- paste(brienz$CHROM, brienz$BIN_START)
neuchatel$loc <- paste(neuchatel$CHROM, neuchatel$BIN_START)

for (i in 1:length(filled$walen_fst)){
  locus <- as.character(filled$loc[i])
  
  info_W <- subset(walen, as.character(walen$loc) == locus)
  filled$walen_fst[i] <- info_W$fst_mean
  
  info_L <- subset(luzern, as.character(luzern$loc) == locus)
  filled$luzern_fst[i] <- info_L$fst_mean
  
  info_B <- subset(brienz, as.character(brienz$loc) == locus)
  filled$brienz_fst[i] <- info_B$fst_mean
  
  info_N <- subset(neuchatel, as.character(neuchatel$loc) == locus)
  filled$neuchatel_fst[i] <- info_N$fst_mean
  
}







#
################# 6.2.1 GO enrichment of lake specific outlier windows ############################### 
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

# GO analysis
library(topGO)

# set node size
nodeS <- 10

# .tsv with 6 columns: "locus.name","relationship","GO.term","GO.ID","aspect","evidence.code"
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
GOraw <- read.csv("Parallel/GO_filt.out", sep = " ", stringsAsFactors = F)

#Lucerne
#outliers <- read.table("./gene_names.txt",header=F,col.names = c("gene"))
outliers <- read.table("Parallel/Parallel_fst/GO_outliers/Lucerne_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))

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
write.table(res.CC_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/L_GOenrichment_CC_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.BP_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/L_GOenrichment_BP_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.MF_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/L_GOenrichment_MF_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)

CC_output <- res.CC_raw_uncorrected
CC_output$class <- "CC"

BP_output <- res.BP_raw_uncorrected
BP_output$class <- "BP"

MF_output <- res.MF_raw_uncorrected
MF_output$class <- "MF"

full_GO_output <- rbind(CC_output, BP_output, MF_output)
write.table(full_GO_output,file="Parallel/Parallel_fst/GO_outliers/L_full_GO_output_nodes10.csv",sep=",",col.names = T,row.names = F)


#Brienz
#outliers <- read.table("./gene_names.txt",header=F,col.names = c("gene"))
outliers <- read.table("Parallel/Parallel_fst/GO_outliers/Brienz_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))

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
write.table(res.CC_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/B_GOenrichment_CC_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.BP_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/B_GOenrichment_BP_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.MF_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/B_GOenrichment_MF_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)

CC_output <- res.CC_raw_uncorrected
CC_output$class <- "CC"

BP_output <- res.BP_raw_uncorrected
BP_output$class <- "BP"

MF_output <- res.MF_raw_uncorrected
MF_output$class <- "MF"

full_GO_output <- rbind(CC_output, BP_output, MF_output)
write.table(full_GO_output,file="Parallel/Parallel_fst/GO_outliers/B_full_GO_output_nodes10.csv",sep=",",col.names = T,row.names = F)

#walen
#outliers <- read.table("./gene_names.txt",header=F,col.names = c("gene"))
outliers <- read.table("Parallel/Parallel_fst/GO_outliers/Walen_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))

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
write.table(res.CC_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/W_GOenrichment_CC_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.BP_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/W_GOenrichment_BP_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.MF_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/W_GOenrichment_MF_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)

CC_output <- res.CC_raw_uncorrected
CC_output$class <- "CC"

BP_output <- res.BP_raw_uncorrected
BP_output$class <- "BP"

MF_output <- res.MF_raw_uncorrected
MF_output$class <- "MF"

full_GO_output <- rbind(CC_output, BP_output, MF_output)
write.table(full_GO_output,file="Parallel/Parallel_fst/GO_outliers/W_full_GO_output_nodes10.csv",sep=",",col.names = T,row.names = F)

#neuchatel
#outliers <- read.table("./gene_names.txt",header=F,col.names = c("gene"))
outliers <- read.table("Parallel/Parallel_fst/GO_outliers/Neuchatel_extractedannotation_unique_gene_names.txt",header=F,col.names = c("gene"))

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
write.table(res.CC_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/N_GOenrichment_CC_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.BP_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/N_GOenrichment_BP_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)
write.table(res.MF_raw,file=paste0("Parallel/Parallel_fst/GO_outliers/N_GOenrichment_MF_nodeSize",nodeS,".csv",sep=""),sep=",",col.names = T,row.names = F)

CC_output <- res.CC_raw_uncorrected
CC_output$class <- "CC"

BP_output <- res.BP_raw_uncorrected
BP_output$class <- "BP"

MF_output <- res.MF_raw_uncorrected
MF_output$class <- "MF"

full_GO_output <- rbind(CC_output, BP_output, MF_output)
write.table(full_GO_output,file="Parallel/Parallel_fst/GO_outliers/N_full_GO_output_nodes10.csv",sep=",",col.names = T,row.names = F)


##now load outliers
luzern_enriched <- read.csv("Parallel/Parallel_fst/GO_outliers/L_full_GO_output_nodes10.csv")
luzern_enriched$lake <- "luzern"
brienz_enriched <- read.csv("Parallel/Parallel_fst/GO_outliers/B_full_GO_output_nodes10.csv")
brienz_enriched$lake <- "brienz"
walen_enriched <- read.csv("Parallel/Parallel_fst/GO_outliers/W_full_GO_output_nodes10.csv")
walen_enriched$lake <- "walen"
neuchatel_enriched <- read.csv("Parallel/Parallel_fst/GO_outliers/N_full_GO_output_nodes10.csv")
neuchatel_enriched$lake <- "neuchatel"

all_enriched <- rbind(luzern_enriched, brienz_enriched, walen_enriched, neuchatel_enriched)

all_enriched$shared <- "no"
for (i in 1:nrow(all_enriched)){
  go_term_name <- as.character(all_enriched$Term)[i]
  go_term_name_df <- subset(all_enriched, as.character(all_enriched$Term) == go_term_name)
  shared_by <- nrow(go_term_name_df)
  all_enriched$shared[i] <- shared_by
}

all_enriched_unique <- all_enriched[!duplicated(all_enriched$Term),]
write.table(all_enriched_unique,file="Parallel/Parallel_fst/GO_outliers/all_shared_enriched_unique.csv",sep=",",col.names = T,row.names = F)
all_enriched_unique_shared <- subset(all_enriched_unique, as.numeric(all_enriched_unique$shared) > 1)
write.table(all_enriched_unique_shared,file="Parallel/Parallel_fst/GO_outliers/all_shared_enriched_unique_shared.csv",sep=",",col.names = T,row.names = F)

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
write.table(all_genes_unique,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique.csv",sep=",",col.names = T,row.names = F)
all_genes_unique_shared <- subset(all_genes_unique, as.numeric(all_genes_unique$shared) > 1)
write.table(all_genes_unique_shared,file="Parallel/Parallel_fst/GO_outliers/all_genes_unique_shared.csv",sep=",",col.names = T,row.names = F)

