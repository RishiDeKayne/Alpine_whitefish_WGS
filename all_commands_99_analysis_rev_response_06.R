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
#write.csv(OUT_brienz, "Parallel/Parallel_fst/Brienz_outliers.csv")
#write.csv(OUT_walen, "Parallel/Parallel_fst/Walen_outliers.csv")
#write.csv(OUT_luzern, "Parallel/Parallel_fst/Lucerne_outliers.csv")
#write.csv(OUT_neuchatel, "Parallel/Parallel_fst/Neuchatel_outliers.csv")

##NEW
min(OUT_brienz$MEAN_FST)
#0.662464
min(OUT_luzern$MEAN_FST)
#0.604515
min(OUT_walen$MEAN_FST)
#0.528866
min(OUT_neuchatel$MEAN_FST)
#0.456683

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
pooled <- read.csv("Parallel/HGR.indivlist_LGR.indivlist_fst.windowed.weir.fst", sep = "\t")
css<- read.csv("Parallel/all_CSS_output_500kpermutations_50000basepair50000step.window.pca.full_noheaders.txt", sep = "\t", head = FALSE)
css2<- read.csv("../../../Downloads/all_CSS_output_100kpermutations_50000basepair50000step.window.pca.full_noheaders.txt", sep = "\t", head = FALSE)
mds<- read.csv("../../../Downloads/all_CSS_output_100kpermutations_50000basepair50000step.window.mds.full_noheaders.txt", sep = "\t", head = FALSE)

mean(css$V5)
mean(css2$V5)
mean(mds$V5)
mean(css$V6)
mean(css2$V6)
mean(mds$V6)

css_colnames6 <- colnames(pooled)
css_colnames5 <- css_colnames6[-c(5)]
colnames(css) <- c(css_colnames5, "q_value")
colnames(css2) <- c(css_colnames5, "q_value")
colnames(mds) <- c(css_colnames5, "q_value")

#now plot css scores vs q value
par(mfrow=c(1,3))
plot(css$MEAN_FST, css$q_value)
plot(css2$MEAN_FST, css2$q_value)
plot(mds$MEAN_FST, mds$q_value)

#css
css_23 <- css[css$N_VARIANTS>23,]
FnE_css<-ecdf(css_23$MEAN_FST)
css_OUT<-cbind(css_23,FST_p=(1-FnE_css(css_23$MEAN_FST)))
OUT_css<-css_OUT[css_OUT$FST_p<0.01,]
OUT_css

css_sig_perm <- subset(css_23, css_23$q_value>0.99)

plot(css_23$MEAN_FST, css_23$q_value, ylab = 'permuted q-value', xlab = "CS'")
points(css_sig_perm$MEAN_FST, css_sig_perm$q_value, col = "blue")
points(OUT_css$MEAN_FST, OUT_css$q_value, col = "red")

#css2
#for plot:
css2$sig <- 1

for (i in 1:nrow(css2)){
  if (css2$q_value[i] > 0.99){
    css2$sig[i] <- 16
  }
}

css2_23 <- css2[css2$N_VARIANTS>23,]
FnE_css2<-ecdf(css2_23$MEAN_FST)
css2_OUT<-cbind(css2_23,FST_p=(1-FnE_css2(css2_23$MEAN_FST)))
OUT_css2<-css2_OUT[css2_OUT$FST_p<0.01,]
OUT_css2

css2_sig_perm <- subset(css2_23, css2_23$q_value>0.99)

plot(css2_23$MEAN_FST, css2_23$q_value, ylab = 'permuted q-value', xlab = "CS'")
points(css2_sig_perm$MEAN_FST, css2_sig_perm$q_value, col = "blue")
points(OUT_css2$MEAN_FST, OUT_css2$q_value, col = "red")


#mds
#for plot:
mds$sig <- 1

for (i in 1:nrow(mds)){
  if (mds$q_value[i] > 0.99){
    mds$sig[i] <- 16
  }
}

mds_23 <- mds[mds$N_VARIANTS>23,]
FnE_mds<-ecdf(mds_23$MEAN_FST)
mds_OUT<-cbind(mds_23,FST_p=(1-FnE_mds(mds_23$MEAN_FST)))
OUT_mds<-mds_OUT[mds_OUT$FST_p<0.01,]
OUT_mds

mds_sig_perm <- subset(mds_23, mds_23$q_value>0.99)

plot(mds_23$MEAN_FST, mds_23$q_value, ylab = 'permuted q-value', xlab = "CS'")
points(mds_sig_perm$MEAN_FST, mds_sig_perm$q_value, col = "blue")
points(OUT_mds$MEAN_FST, OUT_mds$q_value, col = "red")

#and now filter top one percent by css score
#CSS_permuted_outliers <- subset(css, (1-(css$q_value)) < 0.005)
CSS_permuted_outliers <- subset(OUT_css, (1-(OUT_css$q_value)) < 0.01)
CSS2_permuted_outliers <- subset(OUT_css2, (1-(OUT_css2$q_value)) < 0.01)
MDS_permuted_outliers <- subset(OUT_mds, (1-(OUT_mds$q_value)) < 0.01)

par(mfrow=c(3,2))
boxplot(css$MEAN_FST, css2$MEAN_FST, ylab = 'CS score')
axis(1, at = c(1,2), labels = c("old", "new"))
boxplot(css$q_value, css2$q_value, ylab = 'q-value')
axis(1, at = c(1,2), labels = c("old", "new"))

plot(css$MEAN_FST, css2$MEAN_FST, ylab = 'new CS score', xlab = 'old CS score')
plot(css$q_value, css2$q_value, ylab = 'new q-value', xlab = 'old q-value')

plot(css_23$MEAN_FST, css_23$q_value, ylab = 'permuted q-value', xlab = "old CS score")
points(css_sig_perm$MEAN_FST, css_sig_perm$q_value, col = "blue")
points(OUT_css$MEAN_FST, OUT_css$q_value, col = "red")
plot(css2_23$MEAN_FST, css2_23$q_value, ylab = 'permuted q-value', xlab = "new CS score")
points(css2_sig_perm$MEAN_FST, css2_sig_perm$q_value, col = "blue")
points(OUT_css2$MEAN_FST, OUT_css2$q_value, col = "red")

par(mfrow=c(1,1))
plot(css2_23$MEAN_FST, css2_23$q_value, ylab = 'permuted q-value', xlab = "new CS score")
points(css2_sig_perm$MEAN_FST, css2_sig_perm$q_value, col = "blue")
points(OUT_css2$MEAN_FST, OUT_css2$q_value, col = "red")

####see if 1% cs outlier windows are the same for both css2 and mds
test1 <- as.data.frame(cbind(OUT_css2$CHROM, OUT_css2$BIN_END, OUT_css2$MEAN_FST, OUT_css2$q_value))
test1$ID <- paste(test1$V1, test1$V2, sep = "_")
test2 <- as.data.frame(cbind(OUT_mds$CHROM, OUT_mds$BIN_END, OUT_mds$MEAN_FST, OUT_mds$q_value))
test2$ID <- paste(test2$V1, test2$V2, sep = "_")

test3 <- full_join(test1, test2, by = c("ID"))
test4 <- na.omit(test3)

par(mfrow=c(1,2))
plot(css2_23$MEAN_FST, css2_23$BIN_END)
points(css2_sig_perm$MEAN_FST, css2_sig_perm$BIN_END, col = "blue")
points(OUT_css2$MEAN_FST, OUT_css2$BIN_END, col = "red")
points(test4$V3.x, test4$V2.x, col = "green")

plot(mds_23$MEAN_FST, mds_23$BIN_END)
points(mds_sig_perm$MEAN_FST, mds_sig_perm$BIN_END, col = "blue")
points(OUT_mds$MEAN_FST, OUT_mds$BIN_END, col = "red")
points(test4$V3.y, test4$V2.y, col = "green")

####FDR 
#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("qvalue")
library(qvalue)
FDRtest_css2 <- css2_23
FDRtest_css2$p <- 1-FDRtest_css2$q_value
qobj_fdrlevel <- qvalue(p = FDRtest_css2$p, fdr.level = 0.01)
qobj_fdrlevel$significant
FDRtest_css2$sig <- qobj_fdrlevel$significant
true_FDR_css2 <- subset(FDRtest_css2, FDRtest_css2$sig==TRUE)

median(OUT_css2$MEAN_FST)
median(css2_sig_perm$MEAN_FST)
median(true_FDR_css2$MEAN_FST)

FDRtest_mds<- mds_23
FDRtest_mds$p <- 1-FDRtest_mds$q_value
qobj_fdrlevel <- qvalue(p = FDRtest_mds$p, fdr.level = 0.01)
qobj_fdrlevel$significant
FDRtest_mds$sig <- qobj_fdrlevel$significant
true_FDR_mds <- subset(FDRtest_mds, FDRtest_mds$sig==TRUE)

median(OUT_mds$MEAN_FST)
median(mds_sig_perm$MEAN_FST)
median(true_FDR_mds$MEAN_FST)

####see if 1% FDR outlier windows are the same for both css2 and mds
test11 <- as.data.frame(cbind(true_FDR_css2$CHROM, true_FDR_css2$BIN_END, true_FDR_css2$MEAN_FST, true_FDR_css2$q_value))
test11$ID <- paste(test11$V1, test11$V2, sep = "_")
test22 <- as.data.frame(cbind(true_FDR_mds$CHROM, true_FDR_mds$BIN_END, true_FDR_mds$MEAN_FST, true_FDR_mds$q_value))
test22$ID <- paste(test22$V1, test22$V2, sep = "_")

test33 <- full_join(test11, test22, by = c("ID"))
test44 <- na.omit(test33)

par(mfrow=c(1,2))
plot(css2_23$MEAN_FST, css2_23$BIN_END)
points(css2_sig_perm$MEAN_FST, css2_sig_perm$BIN_END, col = "blue")
points(OUT_css2$MEAN_FST, OUT_css2$BIN_END, col = "red")
points(test44$V3.x, test44$V2.x, col = "green")

plot(mds_23$MEAN_FST, mds_23$BIN_END)
points(mds_sig_perm$MEAN_FST, mds_sig_perm$BIN_END, col = "blue")
points(OUT_mds$MEAN_FST, OUT_mds$BIN_END, col = "red")
points(test44$V3.y, test44$V2.y, col = "green")

#write csv of css outliers
#OUTPUT
write.csv(CSS_permuted_outliers,"Parallel/CSS_1percent_outliers.csv", row.names = FALSE)
write.csv(CSS2_permuted_outliers,"Parallel/CSS_1percent_outliers_342_permuted.csv", row.names = FALSE)
write.csv(true_FDR_css2,"Parallel/CSS_1percent_outliers_1659_permuted_FDR.csv", row.names = FALSE)
write.csv(true_FDR_mds,"Parallel/CSS_1percent_outliers_1398_mds_permuted_FDR.csv", row.names = FALSE)


windowsize <- 50000

b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

#pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj.pdf",width=22, height=12)
pdf(file="pca_vs_mds_plot.pdf",width=22, height=12)
par(mfrow=c(2,1),mar=c(0.5,0.15,0,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../99_reanalysis/Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
layout(matrix(seq(1:68),nrow=2,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("slategray1", "slategray2")
altcols <- rep(altcols, as.integer(length(a)/2))

#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css_23$BIN_END[css_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css_23$MEAN_FST[css_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=16,cex=1,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.3),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_css$BIN_END[OUT_css$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css$MEAN_FST[OUT_css$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=1.5)
  #outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  #outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  #points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}

#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css2_23$BIN_END[css2_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css2_23$MEAN_FST[css2_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=as.integer(css2_23$sig),cex=1,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.5),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_css2$BIN_END[OUT_css2$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css2$MEAN_FST[OUT_css2$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=1.5)
  #outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  #outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  #points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
dev.off()


#############Now comparing 342 windows with 1659 from FDR correction
windowsize <- 50000

b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

#pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj.pdf",width=22, height=12)
pdf(file="342_vs_fdr_1659.pdf",width=22, height=12)
par(mfrow=c(2,1),mar=c(0.5,0.15,0,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../99_reanalysis/Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
layout(matrix(seq(1:68),nrow=2,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("slategray1", "slategray2")
altcols <- rep(altcols, as.integer(length(a)/2))

#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css2_23$BIN_END[css2_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css2_23$MEAN_FST[css2_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=as.integer(css2_23$sig),cex=1,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.5),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- OUT_css2$BIN_END[OUT_css2$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css2$MEAN_FST[OUT_css2$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=1.5)
  #outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  #outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  #points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
#fdr
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css2_23$BIN_END[css2_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css2_23$MEAN_FST[css2_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #then plot
  plot(xFst,yFst,pch=as.integer(css2_23$sig),cex=1,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.5),col=altcols[i], cex.axis = 0.5)
  outpoints_xFst <- true_FDR$BIN_END[true_FDR$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- true_FDR$MEAN_FST[true_FDR$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=1.5)
  #outpointsadj_xFst <- test_adj_filt$BIN_END[test_adj_filt$CHROM==a[i]]-(windowsize/2)
  #outpointsadj_yFst <- test_adj_filt$MEAN_FST[test_adj_filt$CHROM==a[i]]
  #points(outpointsadj_xFst,outpointsadj_yFst,col='black',pch=16,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
dev.off()

#############NOW FOR FIG:
windowsize <- 50000

b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

#pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj.pdf",width=22, height=12)
pdf(file="CS_plot_with_significance.pdf",width=15, height=6)
par(mfrow=c(1,1),mar=c(0.5,0.15,0,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../99_reanalysis/Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
layout(matrix(seq(1:34),nrow=1,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols_pooled <- c("grey80", "grey60")
altcols_pooled <- rep(altcols_pooled, as.integer(length(a)/2))
#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css2_23$BIN_END[css2_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css2_23$MEAN_FST[css2_23$CHROM==a[i]]
  #for fst convert values < 0 to 0
  yFst[yFst < 0] <- 0
  #symbol for plot
  #symb <- css2_23$sig[css2_23$CHROM==a[i]]
  #then plot
  plot(xFst,yFst,pch=16,cex=1,ylab="Fst - css",xaxt="n",axes=F,ylim=c(0,0.25),col=altcols_pooled[i], cex.axis = 0.5)
  outpoints_xFst <- css2_sig_perm$BIN_END[css2_sig_perm$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- css2_sig_perm$MEAN_FST[css2_sig_perm$CHROM==a[i]]
  outpoints_yFst[outpoints_yFst < 0] <- 0
  points(outpoints_xFst,outpoints_yFst,col='blue',pch=16,cex=1)
  
  outpoints_xFst <- OUT_css2$BIN_END[OUT_css2$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- OUT_css2$MEAN_FST[OUT_css2$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='red',pch=20,cex=1.5)
  if (i<2){
    axis(side=2)
  }
  #axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
dev.off()

plot(css2_23$MEAN_FST, css2_23$q_value)
points(css2_sig_perm$MEAN_FST, css2_sig_perm$q_value, col = "red")


#############NOW FOR FIG:
windowsize <- 50000

b <- as.vector(chromlist)
#remove chromosome 22, 28, 32, 35, 38
a <- b[-c(22, 28, 32, 35, 38, 40)]

#pdf(file="Parallel//FstPlot_BE_LU_WA_NE_css_permutations_0.01_min10SNPs_3adj.pdf",width=22, height=12)
pdf(file="CS_plot_with_1659_FRD_outliers.pdf",width=12, height=3)
par(mfrow=c(1,1),mar=c(0.5,0.05,0.05,0.15),oma=c(5,4,4,2))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../99_reanalysis/Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
layout(matrix(seq(1:34),nrow=1,ncol=34,byrow=TRUE),width=c(lnew2))

#make colour vector to plot alternating chromosomes alternating colours
altcols_pooled <- c("grey80", "grey60")
altcols_pooled <- rep(altcols_pooled, as.integer(length(a)/2))
#css
for (i in 1:length(a)){
  #extract the x/bp position for each chrom
  xFst <- css2_23$BIN_END[css2_23$CHROM==a[i]]-(windowsize/2)
  #extract the fst values
  yFst <- css2_23$MEAN_FST[css2_23$CHROM==a[i]]
  #then plot
  plot(xFst,yFst,pch=16,cex=0.7,ylab="Fst - css",xaxt="n",axes=F,ylim=c(-0.03,0.25),col=altcols_pooled[i], cex.axis = 0.5)
  outpoints_xFst <- true_FDR_css2$BIN_END[true_FDR_css2$CHROM==a[i]]-(windowsize/2)
  outpoints_yFst <- true_FDR_css2$MEAN_FST[true_FDR_css2$CHROM==a[i]]
  points(outpoints_xFst,outpoints_yFst,col='black',pch=20,cex=0.7)
  if (i<2){
    axis(side=2)
  }
  #axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
}
dev.off()

