###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses
#	1. PCAs
# - 1.1 Full 99 dataset (without outgroups)
# - 1.2 Individual lakes


#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)


################# 1. 1 Full PCA ############################### 
# plot the full 99 dataset (without outgroups)

#load libraries
library(ggplot2)
library("ggrepel")
library(tidyverse)
library(ggplot2)

pca <- read_table2("PCAs/all99_filt_noout.eigenvec", col_names = FALSE)
eigenval <- scan("PCAs/all99_filt_noout.eigenval")

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
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

par(mfrow=c(1,1))

#
tiff("PCAs/plots/99PC1PC2_fig1.tiff", height=8, width=8, units="in", res=300, compression="lzw")

#tiff("PCAs/plots/99PC1PC2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
#tiff("PCAs/plots/99PC1PC2_names.tiff", height=6, width=6, units="in", res=300, compression="lzw")
plot(pca$PC1, 
     pca$PC2, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1))),
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
     ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
     main = "all_samples - P1/P2")
#text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
#       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
#       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
#       pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
dev.off()

tiff("PCAs/plots/99PC1PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
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
#legend('topright', 
       #legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       #col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       #pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
dev.off()

tiff("PCAs/plots/99PC2PC3.tiff", height=8, width=8, units="in", res=300, compression="lzw")
plot(pca$PC2, 
     pca$PC3, 
     col = as.character(pca$colour_plot), 
     pch=pca$symb,
     lwd = 2,
     xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), 
     ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
     main = "all_samples - P2/P3")
#text(pca$PC2, pca$PC3, pca$ind, cex = 0.5, pos = 3)
#legend('topright', 
 #      legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
  #     col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
   #    pt.cex = 1, cex = 0.5, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
dev.off()


###and now to check gill raker/sl associations (confirm findings of CS windows)
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
#Multiple R-squared:  0.06101,	Adjusted R-squared:  0.05034 
#F-statistic: 5.717 on 1 and 88 DF,  p-value: 0.01893

GRC_r2 <- summary(PC1_GRC)$r.squared
GRC_r2_2dp <- format(round(GRC_r2, 3), nsmall = 3)


plot(gill$PC1, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = paste("PC1 vs. gill raker count - R2 =", GRC_r2_2dp))
abline(lm(gill$gill_plot ~ gill$PC1), lty = 3)

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

#and PC1
PC2_GRC <- lm(gill$gill_plot ~ gill$PC2)
summary(PC2_GRC)
#Multiple R-squared:  0.02974,	Adjusted R-squared:  0.01872 
#F-statistic: 2.698 on 1 and 88 DF,  p-value: 0.1041

GRC_r2_PC2 <- summary(PC2_GRC)$r.squared
GRC_r2_PC22dp <- format(round(GRC_r2_PC2, 3), nsmall = 3)


plot(gill$PC2, gill$gill_plot, col = as.character(gill$colour_plot), 
     pch=gill$symb, 
     lwd = 2,
     ylab = "gill raker count",
     xlab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")),
     main = paste("PC2 vs. gill raker count - R2 =", GRC_r2_PC22dp))
#abline(lm(gill$gill_plot ~ gill$PC1), lty = 3)

#
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
tiff("PCAs/all_indivs_PC1_GRC.tiff", height=8, width=10, units="in", res=300, compression="lzw")
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
legend('topright', 
       legend = c((paste("Overall - R2 =", GRC_r2_2dp)),
                  (paste("Overall (exl. profundus) - R2 =", GRC_r2_2dp_NOPROF)),
                  (paste("Brienz/Thun (exl. profundus; non-significant) - R2 =", GRC_B_r2_2dp)), 
                  (paste("Lucerne - R2 =", GRC_L_r2_2dp)), 
                  (paste("Constance (non-significant) - R2 =", GRC_C_r2_2dp)), 
                  (paste("Neuchatel/Biel (non-significant) - R2 =", GRC_N_r2_2dp)), 
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

PC2_GRC <- lm(gill$gill_plot ~ gill$PC2)
summary(PC2_GRC)
PC3_GRC <- lm(gill$gill_plot ~ gill$PC3)
summary(PC3_GRC)
PC4_GRC <- lm(gill$gill_plot ~ gill$PC4)
summary(PC4_GRC)
PC5_GRC <- lm(gill$gill_plot ~ gill$PC5)
summary(PC5_GRC)
PC6_GRC <- lm(gill$gill_plot ~ gill$PC6)
summary(PC6_GRC)
PC7_GRC <- lm(gill$gill_plot ~ gill$PC7)
summary(PC7_GRC)
PC8_GRC <- lm(gill$gill_plot ~ gill$PC8)
summary(PC8_GRC)
PC9_GRC <- lm(gill$gill_plot ~ gill$PC9)
summary(PC9_GRC)
PC10_GRC <- lm(gill$gill_plot ~ gill$PC10)
summary(PC10_GRC)
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
tiff("PCAs/all_indivs_PC1_SL.tiff", height=8, width=10, units="in", res=300, compression="lzw")
plot(SL$PC1, SL$SL_plot, col = as.character(SL$colour_plot), 
     pch=SL$symb, 
     lwd = 2,
     ylab = "SL",
     xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")),
     main = ("PC1 vs. SL"))
abline(lm(SL$SL_plot ~ SL$PC1), lty = 2, col = "black", lwd = 2)
abline(lm(walen_SL$SL_plot ~ walen_SL$PC1), lty = 3, col = "purple4")
#abline(lm(constance_SL$SL_plot ~ constance_SL$PC1), lty = 3, col = "goldenrod")
legend('topright', 
       legend = c((paste("Overall - R2 =", SL_r2_2dp)),
                  (paste("Brienz/Thun (non-significant) - R2 =", SL_B_r2_2dp)), 
                  (paste("Lucerne (non-significant) - R2 =", SL_L_r2_2dp)), 
                  (paste("Constance (non-significant) - R2 =", SL_C_r2_2dp)), 
                  (paste("Neuchatel/Biel (non-significant) - R2 =", SL_N_r2_2dp)), 
                  (paste("Walen/Zurich - R2 =", SL_W_r2_2dp))),
       col = c("black",
               "chartreuse4", 
               "#ff4489", 
               "goldenrod", 
               "turquoise4", 
               "purple4"),
       pt.cex = 1.2, cex = 0.8, pch = "-")
dev.off()



################# 1.2 Individual lake PCAs ############################### 
#load function to extract info and plot lake pca
plotLakePCA <- function(lake, PC1Legend, PC2legend) {
  
  #get lake name
  Lakename <- as.character(lake)
  pca1pos <- as.character(PC1Legend)
  pca2pos <- as.character(PC2legend)
  
  #make eigenvec and eigenval file names
  eigenvec_file <- paste("PCAs/", lake, "_filt_noout.eigenvec", sep = "")
  eigenval_file <- paste("PCAs/", lake, "_filt_noout.eigenval", sep = "")
  
  #eigenvec_file <- "PCAs/Constance_filt_noout.eigenvec"
  #eigenval_file <- "PCAs/Constance_filt_noout.eigenval"
  
  #now load the data
  pca <- read_table2(as.character(eigenvec_file), col_names = FALSE)
  eigenval <- scan(as.character(eigenval_file))
  
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
  
  #now plot
  #tiff(paste("PCAs/plots/", lake, "_P1P2.tiff"),height=8,width=8,units="in",res=200,compression="lzw")
  plot(pca$PC1, 
       pca$PC2, 
       col = as.character(pca$colour_plot), 
       pch=pca$symb,
       lwd = 2,
       xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC2 (", signif(pve$pve[2], 3), "%)")), cex = 1.75,
       main = paste(as.character(lake)))
  #text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend(pca1pos, 
         #legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
         #col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
         #pt.cex = 1.2, cex = 0.6, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24), 
         #y.intersp=0.7, x.intersp=0.7)
  #dev.off()
  
  #now plot
  #tiff(paste("PCAs/plots/", lake, "_P1P3.tiff"),height=8,width=8,units="in",res=200,compression="lzw")
  plot(pca$PC1, 
       pca$PC3, 
       col = as.character(pca$colour_plot), 
       pch=pca$symb,
       lwd = 2,
       xlim = c((min(pca$PC1)-0.005), (max(pca$PC1)+0.1)),
       xlab = (paste0("PC1 (", signif(pve$pve[1], 3), "%)")), 
       ylab = (paste0("PC3 (", signif(pve$pve[3], 3), "%)")), cex = 1.75,
       main = paste(as.character(lake)))
  #text(pca$PC1, pca$PC2, pca$ind, cex = 0.5, pos = 3)
  #legend(pca2pos, 
         #legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Luzern", "Biel", "Neuchatel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
         #col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
         #pt.cex = 1.2, cex = 0.6, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24), 
         #y.intersp=0.7, x.intersp=0.7)
  #dev.off()

  return(print(lake))
}

#can run one lake for e.g.
plotLakePCA("Biel", "bottom", "bottom")
plotLakePCA("Brienz", "bottomright", "bottomright")
plotLakePCA("Thun", "topright", "bottomright")
plotLakePCA("Neuenburg", "bottomright", "bottomright")
plotLakePCA("Lucerne", "bottomright", "topright")
plotLakePCA("Constance", "bottom", "bottom")
plotLakePCA("Zurich", "bottom", "bottom")
plotLakePCA("Walen", "bottom", "bottom")

plotLakePCA("NeuenburgBiel", "bottom", "bottom")
plotLakePCA("ThunBrienz", "topright", "topright")
plotLakePCA("ZurichWalen", "topright", "topright")
















