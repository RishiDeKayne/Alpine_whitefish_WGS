###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	3.  Admixture

#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

fam <- read.csv("Admixture/all99_filt_noout.fam", sep = " ", header = FALSE)

treeorder <- c(34,32,33,64,63,61,21,23,22,51,53,40,38,134,39,55,
               137,93,19,56,16,58,999,17,59,136,129,
               27,26,29,
               128,123458,127,126,132,123470,14,13,12,122,123,131,121,
               109,106,103,104,113,110,111,48,107,70,108,118,114,119,50,66,68,
               47,86,87,90,78,76,79,95,94,96,99,82,83,81,101,98,100,97,102,
               44,41,43,2,5,3,7,72,6,74,71,10)

treeorder_df <- as.data.frame(treeorder)
treeorder_df$order <- 1:91

fam$order <- 0
for (i in 1:nrow(fam)){
  subset_order <- subset(treeorder_df, as.character(treeorder_df$treeorder) == as.character(fam$V1)[i])
  fam$order[i] <- subset_order$order
}

ad2 <- read.table("Admixture/all99_filt_noout.2.Q")
rownames(ad2) <- fam$V1
ad2$order <- fam$order
ad2_o <- ad2[order(ad2$order),]
ad2_of <- ad2_o[,1:2]

ad3 <- read.table("Admixture/all99_filt_noout.3.Q")
rownames(ad3) <- fam$V1
ad3$order <- fam$order
ad3_o <- ad3[order(ad3$order),]
ad3_of <- ad3_o[,1:3]

ad4 <- read.table("Admixture/all99_filt_noout.4.Q")
rownames(ad4) <- fam$V1
ad4$order <- fam$order
ad4_o <- ad4[order(ad4$order),]
ad4_of <- ad4_o[,1:4]

ad5 <- read.table("Admixture/all99_filt_noout.5.Q")
rownames(ad5) <- fam$V1
ad5$order <- fam$order
ad5_o <- ad5[order(ad5$order),]
ad5_of <- ad5_o[,1:5]

ad6 <- read.table("Admixture/all99_filt_noout.6.Q")
rownames(ad6) <- fam$V1
ad6$order <- fam$order
ad6_o <- ad6[order(ad6$order),]
ad6_of <- ad6_o[,1:6]

ad7 <- read.table("Admixture/all99_filt_noout.7.Q")
rownames(ad7) <- fam$V1
ad7$order <- fam$order
ad7_o <- ad7[order(ad7$order),]
ad7_of <- ad7_o[,1:7]

ad8 <- read.table("Admixture/all99_filt_noout.8.Q")
rownames(ad8) <- fam$V1
ad8$order <- fam$order
ad8_o <- ad8[order(ad8$order),]
ad8_of <- ad8_o[,1:8]

ad9 <- read.table("Admixture/all99_filt_noout.9.Q")
rownames(ad9) <- fam$V1
ad9$order <- fam$order
ad9_o <- ad9[order(ad9$order),]
ad9_of <- ad9_o[,1:9]

ad10 <- read.table("Admixture/all99_filt_noout.10.Q")
rownames(ad10) <- fam$V1
ad10$order <- fam$order
ad10_o <- ad10[order(ad10$order),]
ad10_of <- ad10_o[,1:10]

ad11 <- read.table("Admixture/all99_filt_noout.11.Q")
rownames(ad11) <- fam$V1
ad11$order <- fam$order
ad11_o <- ad11[order(ad11$order),]
ad11_of <- ad11_o[,1:11]

ad12 <- read.table("Admixture/all99_filt_noout.12.Q")
rownames(ad12) <- fam$V1
ad12$order <- fam$order
ad12_o <- ad12[order(ad12$order),]
ad12_of <- ad12_o[,1:12]

ad13 <- read.table("Admixture/all99_filt_noout.13.Q")
rownames(ad13) <- fam$V1
ad13$order <- fam$order
ad13_o <- ad13[order(ad13$order),]
ad13_of <- ad13_o[,1:13]

ad14 <- read.table("Admixture/all99_filt_noout.14.Q")
rownames(ad14) <- fam$V1
ad14$order <- fam$order
ad14_o <- ad14[order(ad14$order),]
ad14_of <- ad14_o[,1:14]

tiff("Admixture/test.tiff", height=14, width=14, units="in", res=300, compression="lzw")
par(mfrow = c(13,1))
par(mar=c(3,3,1,1))
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=2", cex.names = 0.3)
barplot(t(as.matrix(ad3_of)), col = c("mediumpurple1", "paleturquoise2", "chartreuse4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=3", cex.names = 0.3)
barplot(t(as.matrix(ad4_of)), col = c("mediumpurple1", "paleturquoise2", "goldenrod", "chartreuse4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=4", cex.names = 0.3)
barplot(t(as.matrix(ad5_of)), col = c("paleturquoise2","goldenrod", "chartreuse2", "mediumpurple1", "chartreuse4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=5", cex.names = 0.3)
barplot(t(as.matrix(ad6_of)), col = c("purple4","mediumpurple1", "goldenrod", "chartreuse4", "#ff4489","paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=6", cex.names = 0.3)
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4","goldenrod", "lightgoldenrod1", "#ff4489", "mediumpurple1","chartreuse2", "paleturquoise2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=7", cex.names = 0.3)
barplot(t(as.matrix(ad8_of)), col = c("mediumpurple1","chartreuse4", "paleturquoise2", "chartreuse2", "#ff4489","goldenrod", "indianred4", "grey"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=8", cex.names = 0.3)
barplot(t(as.matrix(ad9_of)), col = c("chartreuse4","paleturquoise2", "goldenrod", "peachpuff", "#ff4489","purple4", "grey", "mediumpurple1", "chartreuse2"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=9", cex.names = 0.3)
barplot(t(as.matrix(ad10_of)), col = c("purple4","paleturquoise2", "royalblue", "goldenrod", "#ff4489","mediumpurple1", "chartreuse2", "yellow", "chartreuse4", "ivory4"),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=10", cex.names = 0.3)
barplot(t(as.matrix(ad11_of)), col = rainbow(11),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=11", cex.names = 0.3)
barplot(t(as.matrix(ad12_of)), col = rainbow(12),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=12", cex.names = 0.3)
barplot(t(as.matrix(ad13_of)), col = rainbow(13),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=13", cex.names = 0.3)
barplot(t(as.matrix(ad14_of)), col = rainbow(14),
        xlab="individual #", ylab = 'ancestry', border = NA, main = "k=14", cex.names = 0.3)
dev.off()

tiff("Admixture/test2.tiff", height=14, width=14, units="in", res=300, compression="lzw")
par(mfrow = c(9,1))
par(mar=c(1,1,1,1))
barplot(t(as.matrix(ad10_of)), col = c("purple4","paleturquoise2", "royalblue", "goldenrod", "#ff4489","mediumpurple1", "chartreuse2", "yellow", "chartreuse4", "ivory4"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad9_of)), col = c("chartreuse4","paleturquoise2", "goldenrod", "peachpuff", "#ff4489","purple4", "grey", "mediumpurple1", "chartreuse2"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad8_of)), col = c("mediumpurple1","chartreuse4", "paleturquoise2", "chartreuse2", "#ff4489","goldenrod", "indianred4", "grey"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad7_of)), col = c("chartreuse4","goldenrod", "lightgoldenrod1", "#ff4489", "mediumpurple1","chartreuse2", "paleturquoise2"),
        xlab="individual #", ylab = 'm', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad6_of)), col = c("purple4","mediumpurple1", "goldenrod", "chartreuse4", "#ff4489","paleturquoise2"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad5_of)), col = c("paleturquoise2","goldenrod", "chartreuse4", "mediumpurple1", "chartreuse2"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad4_of)), col = c("mediumpurple1", "paleturquoise2", "goldenrod", "chartreuse4"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad3_of)), col = c("mediumpurple1", "paleturquoise2", "chartreuse4"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
barplot(t(as.matrix(ad2_of)), col = c("chartreuse4", "paleturquoise2"),
        xlab="individual #", ylab = 'n', cex.names = 0.3, xaxt="n", yaxt="n")
dev.off()

#and plot error for each k:
CV_error <- read.csv("Admixture/admixture_cv_error_output.txt", sep = "\t", header = FALSE)

tiff("Admixture/CV_plot.tiff", height=4, width=4, units="in", res=300, compression="lzw")
plot(CV_error$V1, CV_error$V2, pch = 16, ylab = "CV error from admixture run", xlab = "K")
lines(CV_error$V1, CV_error$V2, col = "grey50", lty = 3)
dev.off()
