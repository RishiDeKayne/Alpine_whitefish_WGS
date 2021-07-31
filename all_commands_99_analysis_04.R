###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	4.  Dstats
# 4.1 - Large pelagic introgression?
# 4.2 - Profundal - C. profundus and C. nobilis introgression?
# 4.3 - large pleagic acrinasus - Thun introgression

#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

################# 4. 1 Large Pelagic ############################### 
library(vioplot)
library(ggplot2)
library(gridExtra)
#dsuites for pelagic int
full <- read.csv("Dstats/Luzern_filt.csv", sep = ",", header = F)
colnames(full) <- c("P1", "P2", "P3", "Dstatistic", "p-value", "f_G")

full_nonobilis <- full[-grep("nobilis", full$P1),]
full_nonobilis <- full_nonobilis[-grep("nobilis", full_nonobilis$P2),]

full_no047 <- full_nonobilis[-grep("047", full_nonobilis$P1),]
full_no047 <- full_no047[-grep("047", full_no047$P2),]

#now isolate pelagic individuals
pelagicv1 <- full_no047[grep("pelagic", full_no047$P1),]
pelagicv2 <- full_no047[grep("pelagic", full_no047$P2),]

pelagic <- rbind(pelagicv1, pelagicv2)

pelagic$duplicate <- "no"
for(i in 1:length(pelagic$P1)){
  trio <- pelagic[i,]
  a <- nrow(trio[grep("pelagic", trio$P1),])
  b <- nrow(trio[grep("pelagic", trio$P2),])
  if(as.integer(a) + as.integer(b) > 1){
    pelagic$duplicate[i] <- 'yes'
  }
}

pelagic_nodup <- subset(pelagic, pelagic$duplicate == "no")

##filter out pelagic topologies from full
full_no047_nopel <- full_no047[-grep("pelagic", full_no047$P1),]
full_no047_nopel <- full_no047_nopel[-grep("pelagic", full_no047_nopel$P2),]

full_no047_nopel_dups <- full_no047_nopel
#want to remove double C.nobilis - C.zugensis - C.bodenbalchen - C.benthic_intermediate - C.alpnacherfelchen
full_no047_nopel_dups$duplicate <- "no"
for(i in 1:length(full_no047_nopel_dups$P1)){
  trio <- full_no047_nopel_dups[i,]
  a <- nrow(trio[grep("nobilis", trio$P1),])
  b <- nrow(trio[grep("nobilis", trio$P2),])
  c <- nrow(trio[grep("zugensis", trio$P1),])
  d <- nrow(trio[grep("zugensis", trio$P2),])  
  e <- nrow(trio[grep("bodenbalchen", trio$P1),])
  f <- nrow(trio[grep("bodenbalchen", trio$P2),]) 
  g <- nrow(trio[grep("benthic_intermediate", trio$P1),])
  h <- nrow(trio[grep("benthic_intermediate", trio$P2),])
  j <- nrow(trio[grep("alpnacherfelchen", trio$P1),])
  k <- nrow(trio[grep("alpnacherfelchen", trio$P2),])
  if(as.integer(a) + as.integer(b) > 1){
    full_no047_nopel_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(c) + as.integer(d) > 1){
    full_no047_nopel_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(e) + as.integer(f) > 1){
    full_no047_nopel_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(g) + as.integer(h) > 1){
    full_no047_nopel_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(j) + as.integer(k) > 1){
    full_no047_nopel_dups$duplicate[i] <- 'yes'
  }
}

full_no047_nopel <- subset(full_no047_nopel_dups, as.character(full_no047_nopel_dups$duplicate) == "no") 

#plot full without pelagic vs. pelagic (duplicates i.e. P1 and P2 are both pelagic)
par(mfrow=c(2,1))

boxplot(full_no047_nopel$Dstatistic, pelagic_nodup$Dstatistic,
        names=c("Lucerne (constance)", "pelagic (constance)"),
        cex.axis = 0.7,
        ylab = "Dmin",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("plum3", "royalblue3"))

boxplot(-log10(full_no047_nopel$`p-value`), -log10(pelagic_nodup$`p-value`),
        names=c("Lucerne (constance)", "pelagic (constance)"),
        cex.axis = 0.7,
        ylab = "-log10 Pvalue",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("plum3", "royalblue3"))

wilcox.test(full_no047_nopel$Dstatistic, pelagic_nodup$Dstatistic)
#W = 30145, p-value = 9.034e-14
#alternative hypothesis: true location shift is not equal to 0

#separate full file by wartmanii vs. other constance spp
full_no047_nopel_wartmanii <- full_no047_nopel[grep("wartmanii", full_no047_nopel$P3),]
full_no047_nopel_otherconstance <- full_no047_nopel[-grep("wartmanii", full_no047_nopel$P3),]

full_no047_nopel_macrophthalmus <- full_no047_nopel_otherconstance[grep("C.macrophthalmus", full_no047_nopel_otherconstance$P3),]
full_no047_nopel_arenicolus <- full_no047_nopel_otherconstance[-grep("C.macrophthalmus", full_no047_nopel_otherconstance$P3),]

#separate pelagic no duplicates file by wartmanii vs. other constance spp
pelagic_nodup_wartmanii <- pelagic_nodup[grep("wartmanii", pelagic_nodup$P3),]
pelagic_nodup_otherconstance <- pelagic_nodup[-grep("wartmanii", pelagic_nodup$P3),]

pelagic_nodup_macrophthalmus <- pelagic_nodup_otherconstance[grep("C.macrophthalmus", pelagic_nodup_otherconstance$P3),]
pelagic_nodup_arenicolus <- pelagic_nodup_otherconstance[-grep("C.macrophthalmus", pelagic_nodup_otherconstance$P3),]

par(mfrow=c(2,1))
boxplot(full_no047_nopel_wartmanii$Dstatistic, full_no047_nopel_macrophthalmus$Dstatistic, full_no047_nopel_arenicolus$Dstatistic, 
        pelagic_nodup_wartmanii$Dstatistic, pelagic_nodup_macrophthalmus$Dstatistic, pelagic_nodup_arenicolus$Dstatistic,
        names=c("Lucerne (wart)", "Lucerne (marcroph)", "Lucerne (arenic)",
                "pelagic (wart)", "pelagic (marcroph)","pelagic (arenic)"),
        cex.axis = 0.7,
        ylab = "Dmin",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("tan4", "goldenrod3", "goldenrod1", "tan4", "goldenrod3", "goldenrod1"))

boxplot(-log10(full_no047_nopel_wartmanii$`p-value`), -log10(full_no047_nopel_macrophthalmus$`p-value`), -log10(full_no047_nopel_arenicolus$`p-value`), 
        -log10(pelagic_nodup_wartmanii$`p-value`), -log10(pelagic_nodup_macrophthalmus$`p-value`), -log10(pelagic_nodup_arenicolus$`p-value`),
        names=c("Lucerne (wart)", "Lucerne (marcroph)", "Lucerne (arenic)",
                "pelagic (wart)", "pelagic (marcroph)","pelagic (arenic)"),
        cex.axis = 0.7,
        ylab = "-log10(p-value)",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("tan4", "goldenrod3", "goldenrod1", "tan4", "goldenrod3", "goldenrod1"))
abline(h=2, col = "red", lty = 2)

library(vioplot)
vioplot(full_no047_nopel_wartmanii$Dstatistic, full_no047_nopel_macrophthalmus$Dstatistic, full_no047_nopel_arenicolus$Dstatistic, 
        pelagic_nodup_wartmanii$Dstatistic, pelagic_nodup_macrophthalmus$Dstatistic, pelagic_nodup_arenicolus$Dstatistic,
        names=c("Lucerne (wart)", "Lucerne (marcroph)", "Lucerne (arenic)",
                "pelagic (wart)", "pelagic (marcroph)","pelagic (arenic)"),
        cex.axis = 0.7,
        ylab = "Dmin",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("tan4", "goldenrod3", "goldenrod1", "tan4", "goldenrod3", "goldenrod1"))

vioplot(-log10(full_no047_nopel_wartmanii$`p-value`), -log10(full_no047_nopel_macrophthalmus$`p-value`), -log10(full_no047_nopel_arenicolus$`p-value`), 
        -log10(pelagic_nodup_wartmanii$`p-value`), -log10(pelagic_nodup_macrophthalmus$`p-value`), -log10(pelagic_nodup_arenicolus$`p-value`),
        names=c("Lucerne (wart)", "Lucerne (marcroph)", "Lucerne (arenic)",
                "pelagic (wart)", "pelagic (marcroph)","pelagic (arenic)"),
        cex.axis = 0.7,
        ylab = "-log10(p-value)",
        main = "All Lucerne P1 and P2 vs. pelagic int",
        col = c("tan4", "goldenrod3", "goldenrod1", "tan4", "goldenrod3", "goldenrod1"))
abline(h=2, col = "red", lty = 2)

####################now plot with points
full_no047_nopel_wartmanii$group <- "aFU_w"
full_no047_nopel_macrophthalmus$group <- "cFU_m"
full_no047_nopel_arenicolus$group <- "eFU_a"
pelagic_nodup_wartmanii$group <- "bPE_w"
pelagic_nodup_macrophthalmus$group <- "dPE_m"
pelagic_nodup_arenicolus$group  <- "fPE_a"

test_all <- rbind(full_no047_nopel_wartmanii, full_no047_nopel_macrophthalmus, full_no047_nopel_arenicolus, 
                  pelagic_nodup_wartmanii, pelagic_nodup_macrophthalmus, pelagic_nodup_arenicolus)

test_all$log <- -log10(test_all$`p-value`)

test_all$sig <- "no"
for(i in 1:nrow(test_all)){
  if(test_all$log[i] >= 2){
    test_all$sig[i] <- "yes"
  }
}

exp1 <- expression(paste("All - ", italic("C. wartmanni")))
exp2 <- expression(paste("LP - ", italic("C. wartmanni")))
exp3 <- expression(paste("All - ", italic("C. macrophthalmus")))
exp4 <- expression(paste("LP - ", italic("C. macrophthalmus")))
exp5 <- expression(paste("All - ", italic("C. arenicolus")))
exp6 <- expression(paste("LP - ", italic("C. arenicolus")))

library(ggsignif)

ggplot(data=test_all,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("tan4", "tan4", "goldenrod3", "goldenrod3", "goldenrod1", "goldenrod1"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,16)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("darkgrey", "black", "darkgrey", "black", "darkgrey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0)+
  scale_x_discrete(labels=c("aFU_w" = exp1, "bPE_w" = exp2, "cFU_m" = exp3,
                            "dPE_m" = exp4, "eFU_a" = exp5, "fPE_a" = exp6))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_w", "bPE_w")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  geom_signif(comparisons=list(c("cFU_m", "dPE_m")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()

b <- ggplot(data=test_all,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("tan4", "tan4", "goldenrod3", "goldenrod3", "goldenrod1", "goldenrod1"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,16)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("darkgrey", "black", "darkgrey", "black", "darkgrey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0)+
  scale_x_discrete(labels=c("aFU_w" = exp1, "bPE_w" = exp2, "cFU_m" = exp3,
                            "dPE_m" = exp4, "eFU_a" = exp5, "fPE_a" = exp6))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_w", "bPE_w")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  geom_signif(comparisons=list(c("cFU_m", "dPE_m")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()


#wartmanni test in each
wilcox.test(full_no047_nopel_wartmanii$Dstatistic, pelagic_nodup_wartmanii$Dstatistic)
#W = 224, p-value < 2.2e-16
##SIG DIFF

#macrophthalmus test in each
wilcox.test(full_no047_nopel_macrophthalmus$Dstatistic, pelagic_nodup_macrophthalmus$Dstatistic)
#W = 3281, p-value = 0.0002781
##SIG DIFF

#arenicolus test in each
wilcox.test(full_no047_nopel_arenicolus$Dstatistic, pelagic_nodup_arenicolus$Dstatistic)
#W = 6173, p-value = 0.1663
#NOT SIG DIFF

wilcox.test(full_no047_nopel_wartmanii$Dstatistic, pelagic_nodup_wartmanii$Dstatistic, alternative = "less")
wilcox.test(full_no047_nopel_macrophthalmus$Dstatistic, pelagic_nodup_macrophthalmus$Dstatistic, alternative = "less")
wilcox.test(full_no047_nopel_arenicolus$Dstatistic, pelagic_nodup_arenicolus$Dstatistic, alternative = "less")


################# 4. 2 Profundal ############################### 

## 4.2.1 dsuites for profundal from Thun

#load full dataset of Thun with Biel/Neuchatel as P3
full <- read.csv("Dstats/ThunvsBielNeuchatel.csv", sep = ",", header = F)
colnames(full) <- c("P1", "P2", "P3", "Dstatistic", "p-value", "f_G")

#remove acrinasus since it is very different/doesnt cluster with Thun anyway
full_noacrinasus <- full[-grep("C.acrinasus", full$P1),]
full_noacrinasus <- full_noacrinasus[-grep("C.acrinasus", full_noacrinasus$P2),]

#now isolate topologies where profundus individuals are P1 or P2
profv1 <- full_noacrinasus[grep("C.profundus", full_noacrinasus$P1),]
profv2 <- full_noacrinasus[grep("C.profundus", full_noacrinasus$P2),]

#combine
prof <- rbind(profv1, profv2)

#remove rows where both P1 and P2 are profundus
prof$duplicate <- "no"
for(i in 1:length(prof$P1)){
  trio <- prof[i,]
  a <- nrow(trio[grep("C.profundus", trio$P1),])
  b <- nrow(trio[grep("C.profundus", trio$P2),])
  if(as.integer(a) + as.integer(b) > 1){
    prof$duplicate[i] <- 'yes'
  }
}

prof_nodup_T <- subset(prof, prof$duplicate == "no")

#then remove all profundus rows from full df
full_noacrinasus_noprof <- full_noacrinasus[-grep("C.profundus", full_noacrinasus$P1),]
full_noacrinasus_noprof <- full_noacrinasus_noprof[-grep("C.profundus", full_noacrinasus_noprof$P2),]

full_noacrinasus_noprof_dups <- full_noacrinasus_noprof
#want to remove double C.nobilis - C.zugensis - C.bodenbalchen - C.benthic_intermediate - C.alpnacherfelchen
full_noacrinasus_noprof_dups$duplicate <- "no"
for(i in 1:length(full_noacrinasus_noprof_dups$P1)){
  trio <- full_noacrinasus_noprof_dups[i,]
  c <- nrow(trio[grep("alpinus", trio$P1),])
  d <- nrow(trio[grep("alpinus", trio$P2),])  
  e <- nrow(trio[grep("fatioi", trio$P1),])
  f <- nrow(trio[grep("fatioi", trio$P2),]) 
  g <- nrow(trio[grep("steinmanni", trio$P1),])
  h <- nrow(trio[grep("steinmanni", trio$P2),])
  j <- nrow(trio[grep("albellus", trio$P1),])
  k <- nrow(trio[grep("albellus", trio$P2),])
  if(as.integer(c) + as.integer(d) > 1){
    full_noacrinasus_noprof_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(e) + as.integer(f) > 1){
    full_noacrinasus_noprof_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(g) + as.integer(h) > 1){
    full_noacrinasus_noprof_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(j) + as.integer(k) > 1){
    full_noacrinasus_noprof_dups$duplicate[i] <- 'yes'
  }
}

full_noacrinasus_noprof <- subset(full_noacrinasus_noprof_dups, as.character(full_noacrinasus_noprof_dups$duplicate) == "no") 


#and plot comparison
boxplot(full_noacrinasus_noprof$Dstatistic, prof_nodup_T$Dstatistic,
        names=c("all Thun", "C.profundus"),
        cex.axis = 1,
        ylab = "Dmin",
        main = "All Thun P1 and P2 vs. profundus - P3 = Biel/Neuchatel",
        col = "lightblue")

boxplot(-log10(full_noacrinasus_noprof$`p-value`), -log10(prof_nodup_T$`p-value`),
        names=c("all Thun", "C.profundus"),
        cex.axis = 1,
        ylab = "-log10(p-value)",
        main = "All Thun P1 and P2 vs. profundus - P3 = Biel/Neuchatel")


## 4.2.2 dsuites for nobilis from Luzern
full <- read.csv("Dstats/LuzernvsAlbellus_filt.csv", sep = ",", header = F)
colnames(full) <- c("P1", "P2", "P3", "Dstatistic", "p-value", "f_G")

#then remove pelagic int we think is a pelagic int from df
full_nopelagic <- full[-grep("C.pelagic_intermediate", full$P1),]
full_nopelagic <- full_nopelagic[-grep("C.pelagic_intermediate", full_nopelagic$P2),]

full_no047 <- full_nopelagic[-grep("047", full_nopelagic$P1),]
full_no047 <- full_no047[-grep("047", full_no047$P2),]

#now isolate topologies with nobilis as P1 or P2
profv1 <- full_no047[grep("C.nobilis", full_no047$P1),]
profv2 <- full_no047[grep("C.nobilis", full_no047$P2),]

prof <- rbind(profv1, profv2)

#again remove topologies where both P1 and P2 are nobilis
prof$duplicate <- "no"
for(i in 1:length(prof$P1)){
  trio <- prof[i,]
  a <- nrow(trio[grep("C.nobilis", trio$P1),])
  b <- nrow(trio[grep("C.nobilis", trio$P2),])
  if(as.integer(a) + as.integer(b) > 1){
    prof$duplicate[i] <- 'yes'
  }
}

prof_nodup_L <- subset(prof, prof$duplicate == "no")

#and remove nobilis from full df
full_nopelagic_nonobilis <- full_no047[-grep("C.nobilis", full_no047$P1),]
full_nopelagic_nonobilis <- full_nopelagic_nonobilis[-grep("C.nobilis", full_nopelagic_nonobilis$P2),]

full_nopelagic_nonobilis_dups <- full_nopelagic_nonobilis
#want to remove double C.nobilis - C.zugensis - C.bodenbalchen - C.benthic_intermediate - C.alpnacherfelchen
full_nopelagic_nonobilis_dups$duplicate <- "no"
for(i in 1:length(full_nopelagic_nonobilis_dups$P1)){
  trio <- full_nopelagic_nonobilis_dups[i,]
  c <- nrow(trio[grep("zugensis", trio$P1),])
  d <- nrow(trio[grep("zugensis", trio$P2),])  
  e <- nrow(trio[grep("bodenbalchen", trio$P1),])
  f <- nrow(trio[grep("bodenbalchen", trio$P2),]) 
  g <- nrow(trio[grep("benthic_intermediate", trio$P1),])
  h <- nrow(trio[grep("benthic_intermediate", trio$P2),])
  j <- nrow(trio[grep("alpnacherfelchen", trio$P1),])
  k <- nrow(trio[grep("alpnacherfelchen", trio$P2),])
  if(as.integer(c) + as.integer(d) > 1){
    full_nopelagic_nonobilis_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(e) + as.integer(f) > 1){
    full_nopelagic_nonobilis_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(g) + as.integer(h) > 1){
    full_nopelagic_nonobilis_dups$duplicate[i] <- 'yes'
  }
  if(as.integer(j) + as.integer(k) > 1){
    full_nopelagic_nonobilis_dups$duplicate[i] <- 'yes'
  }
}

full_nopelagic_nonobilis <- subset(full_nopelagic_nonobilis_dups, as.character(full_nopelagic_nonobilis_dups$duplicate) == "no") 


#plot luzern
boxplot(full_nopelagic_nonobilis$Dstatistic, prof_nodup_L$Dstatistic,
        names=c("all Luzern", "C.nobilis"),
        cex.axis = 1,
        ylab = "Dmin",
        main = "All Luzern P1 and P2 vs. nobilis - P3 = C.albellus",
        col = "darkgreen")

boxplot(-log10(full_nopelagic_nonobilis$`p-value`), -log10(prof_nodup_L$`p-value`),
        names=c("all Luzern", "C.nobilis"),
        cex.axis = 1,
        ylab = "-log10(p-value)",
        main = "All Luzern P1 and P2 vs. nobilis - P3 = C.albellus")


######################################
#now plot all together showing both profundal morphs are different
#C.nobilis has introgression from Thun/Brienz (C. albellus) and profundus seems to have no excess alleles with Biel/Neuchatel
boxplot(full_nopelagic_nonobilis$Dstatistic, prof_nodup_L$Dstatistic, full_noacrinasus_noprof$Dstatistic, prof_nodup_T$Dstatistic,
        names=c("all Luzern", "C.nobilis", "all Thun", "C.profundus"),
        cex.axis = 1,
        ylab = "Dmin",
        col = c("darkgreen", "darkgreen", "lightblue", "lightblue"))

boxplot(-log10(full_nopelagic_nonobilis$`p-value`), -log10(prof_nodup_L$`p-value`), -log10(full_noacrinasus_noprof$`p-value`), -log10(prof_nodup_T$`p-value`),
        names=c("all Luzern", "C.nobilis", "all Thun", "C.profundus"),
        cex.axis = 1,
        ylab = "-log10(p-value)",
        col = c("darkgreen", "darkgreen", "lightblue", "lightblue"))
abline(h=2, col = "red")

wilcox.test(full_nopelagic_nonobilis$Dstatistic, prof_nodup_L$Dstatistic)
#W = 70, p-value = 1.422e-11
wilcox.test(full_noacrinasus_noprof$Dstatistic, prof_nodup_T$Dstatistic)
#W = 247032, p-value = 0.0001335

####################now plot with points
full_nopelagic_nonobilis$group <- "aFU_L"
prof_nodup_L$group <- "bPR_L"
full_noacrinasus_noprof$group <- "cFU_T"
prof_nodup_T$group <- "dPR_T"


test_all2 <- rbind(full_nopelagic_nonobilis, prof_nodup_L, full_noacrinasus_noprof, prof_nodup_T)

test_all2$log <- -log10(test_all2$`p-value`)

test_all2$sig <- "no"
for(i in 1:nrow(test_all2)){
  if(test_all2$log[i] >= 2){
    test_all2$sig[i] <- "yes"
  }
  if5070 <- test_all2[i,]
  e <- nrow(if5070[grep("070_C.alpnacherfelchen_Lucerne_26004", if5070$P1),])
  f <- nrow(if5070[grep("070_C.alpnacherfelchen_Lucerne_26004", if5070$P2),])
  if(as.integer(e) + as.integer(f) > 0 && test_all2$log[i] >= 2){
    test_all2$sig[i] <- "yes70"
  }
  if(as.integer(e) + as.integer(f) > 0 && test_all2$log[i] < 2){
    test_all2$sig[i] <- "no70"
  }
}

exp7 <- expression(paste("All Luzern - ", italic("C.albellus")))
exp8 <- expression(paste(italic("C. nobilis"), " - ", italic("C.albellus")))
exp9 <- expression(paste("All Thun - Biel/Neuchâtel"))
exp10 <- expression(paste(italic("C. profundus"), " - Biel/Neuchâtel"))


#quartz()
library(ggplot2)
ggplot(data=test_all2,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("chartreuse4", "chartreuse4", "turquoise4", "turquoise4"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,2,16,17)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("grey", "black", "grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0)+
  scale_x_discrete(labels=c("aFU_L" = exp7, "bPR_L" = exp8, 
                            "cFU_T" = exp9, "dPR_T" = exp10))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_L", "bPR_L")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()


test_all_luzern <- rbind(full_nopelagic_nonobilis, prof_nodup_L)

test_all_luzern$log <- -log10(test_all_luzern$`p-value`)

test_all_luzern$sig <- "no"
for(i in 1:nrow(test_all_luzern)){
  if(test_all_luzern$log[i] >= 2){
    test_all_luzern$sig[i] <- "yes"
  }
  if5070 <- test_all_luzern[i,]
  e <- nrow(if5070[grep("070_C.alpnacherfelchen_Lucerne_26004", if5070$P1),])
  f <- nrow(if5070[grep("070_C.alpnacherfelchen_Lucerne_26004", if5070$P2),])
  if(as.integer(e) + as.integer(f) > 0 && test_all_luzern$log[i] >= 2){
    test_all_luzern$sig[i] <- "yes70"
  }
  if(as.integer(e) + as.integer(f) > 0 && test_all_luzern$log[i] < 2){
    test_all_luzern$sig[i] <- "no70"
  }
}

a <- ggplot(data=test_all_luzern,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("chartreuse4", "chartreuse4"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,2,16,17)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("grey", "black", "grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0)+
  scale_x_discrete(labels=c("aFU_L" = exp7, "bPR_L" = exp8))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_L", "bPR_L")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()
  theme(legend.position = "none")


test_all_thun <- rbind(full_noacrinasus_noprof, prof_nodup_T)

test_all_thun$log <- -log10(test_all_thun$`p-value`)

test_all_thun$sig <- "no"
for(i in 1:nrow(test_all_thun)){
  if(test_all_thun$log[i] >= 2){
    test_all_thun$sig[i] <- "yes"
  }
}

c <- ggplot(data=test_all_thun,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("turquoise4", "turquoise4"),lwd=2, alpha = 0.4, show.legend = FALSE)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,16)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0)+
  scale_x_discrete(labels=c("cFU_T" = exp9, "dPR_T" = exp10))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("cFU_T", "dPR_T")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()
  theme(legend.position = "none")

  
#####
aa <- ggplot(data=test_all_luzern,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("chartreuse4", "chartreuse4"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,2,16,17)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("grey", "black", "grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0, show.legend = FALSE)+
  scale_x_discrete(labels=c("aFU_L" = exp7, "bPR_L" = exp8))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_L", "bPR_L")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()

bb <- ggplot(data=test_all,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("tan4", "tan4", "goldenrod3", "goldenrod3", "goldenrod1", "goldenrod1"),lwd=2, alpha = 0.4,show.legend = F)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,16)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("darkgrey", "black", "darkgrey", "black", "darkgrey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0, show.legend = FALSE)+
  scale_x_discrete(labels=c("aFU_w" = exp1, "bPE_w" = exp2, "cFU_m" = exp3,
                            "dPE_m" = exp4, "eFU_a" = exp5, "fPE_a" = exp6))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  geom_signif(comparisons=list(c("aFU_w", "bPE_w")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  geom_signif(comparisons=list(c("cFU_m", "dPE_m")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()

cc <- ggplot(data=test_all_thun,aes(x=group, y=Dstatistic)) + 
  #geom_point(aes(color=group,shape=group), na.rm=TRUE) +  
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("turquoise4", "turquoise4"),lwd=2, alpha = 0.4, show.legend = FALSE)+
  #geom_violin(alpha=0.8, position = position_dodge(width = .75),size=1,color=NA) +
  scale_shape_manual(values=c(1,16)) + 
  #scale_color_manual(values=c("grey", "black", "grey", "black", "grey", "black"))+
  scale_color_manual(values=c("grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0, show.legend = FALSE)+
  scale_x_discrete(labels=c("cFU_T" = exp9, "dPR_T" = exp10))+
  #geom_hline(yintercept=2, linetype="dashed", color = "red")+
  #geom_signif(comparisons=list(c("cFU_T", "dPR_T")), annotations="***",
  #            y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()

#plots for FIGURE 2
tiff("Figures/Figure2.tiff", height=8, width=8, units="in", res=300, compression="lzw")
grid.arrange(bb,                                      # bar plot spaning two columns
             aa, 
             cc,                               # box plot and scatter plot
             ncol = 1, nrow = 3, 
             layout_matrix = rbind(c(1,1), c(2,3)))
dev.off()

################# 4. 3 C.acrinasus thun ############################### 
#check acrinasus introgression from Thun:
#acr_con_thun_filt_all_noprof.csv and Con_Con_Thun_filt_noprof.csv
full_acr <- read.csv("Dstats/acr_con_thun_filt_all_noprof.csv", sep = ",", header = F)
colnames(full_acr) <- c("P1", "P2", "P3", "Dstatistic", "p-value", "f_G")

#remove topologies where eitehr P1 and P3 or P2 and P3 are acrinasus
full_acr$duplicate1 <- "no"
for(i in 1:length(full_acr$P1)){
  trio <- full_acr[i,]
  a <- nrow(trio[grep("C.acrinasus", trio$P1),])
  b <- nrow(trio[grep("C.acrinasus", trio$P3),])
  if(as.integer(a) + as.integer(b) > 1){
    full_acr$duplicate1[i] <- 'yes'
  }
}

full_acr$duplicate2 <- "no"
for(i in 1:length(full_acr$P1)){
  trio <- full_acr[i,]
  a <- nrow(trio[grep("C.acrinasus", trio$P2),])
  b <- nrow(trio[grep("C.acrinasus", trio$P3),])
  if(as.integer(a) + as.integer(b) > 1){
    full_acr$duplicate2[i] <- 'yes'
  }
}

full_acr_nodups <- subset(full_acr, as.character(full_acr$duplicate1) == "no") 
full_acr_nodups2 <- subset(full_acr_nodups, as.character(full_acr_nodups$duplicate2) == "no") 

full_con <- read.csv("Dstats/Con_Con_Thun_filt_noprof.csv", sep = ",", header = F)
colnames(full_con) <- c("P1", "P2", "P3", "Dstatistic", "p-value", "f_G")

#again remove topologies where both P1 and P2 are either macrophthalmus or both arenicolus
full_con$duplicate1 <- "no"
full_con$duplicate2 <- "no"

for(i in 1:length(full_con$P1)){
  trio <- full_con[i,]
  a <- nrow(trio[grep("C.macrophthalmus", trio$P1),])
  b <- nrow(trio[grep("C.macrophthalmus", trio$P2),])
  c <- nrow(trio[grep("C.arenicolus", trio$P1),])
  d <- nrow(trio[grep("C.arenicolus", trio$P2),])
  e <- nrow(trio[grep("C.wartmanii", trio$P1),])
  f <- nrow(trio[grep("C.wartmanii", trio$P2),])
  if(as.integer(a) + as.integer(b) > 1){
    full_con$duplicate1[i] <- 'yes'
  }
  if(as.integer(c) + as.integer(d) > 1){
    full_con$duplicate1[i] <- 'yes'
  }
  if(as.integer(e) + as.integer(f) > 1){
    full_con$duplicate1[i] <- 'yes'
  }
}

full_con_nodups <- subset(full_con, as.character(full_con$duplicate1) == "no") 

#remove topologies with acrinasus at P3
full_con_nodups <- full_con_nodups[-grep("C.acrinasus", full_con_nodups$P3),]


#plot luzern
boxplot(full_con_nodups$Dstatistic, full_acr_nodups2$Dstatistic,
        names=c("macrophthalmus/arenicolus", "C.acrinasus"),
        cex.axis = 1,
        ylab = "Dmin",
        main = "macrophthalmus, arenicolus P1 and P2 vs. acrinasus - P3 = Thun-profundus",
        col = "darkgreen")

wilcox.test(full_con_nodups$Dstatistic, full_acr_nodups2$Dstatistic)
mean(full_acr_nodups2$Dstatistic)

full_acr_nodups2$group <- "g2"
full_con_nodups$group <- "g1"

all_con_acr <- rbind(full_con_nodups, full_acr_nodups2)

all_con_acr$log <- -log10(all_con_acr$`p-value`)

all_con_acr$sig <- "no"
for(i in 1:nrow(all_con_acr)){
  if(all_con_acr$log[i] >= 2){
    all_con_acr$sig[i] <- "yes"
  }
}

dd <- ggplot(data=all_con_acr,aes(x=group, y=Dstatistic)) + 
  geom_boxplot(notch = TRUE,  outlier.size = -1, color=c("chartreuse4", "chartreuse4"),lwd=2, alpha = 0.4, show.legend = FALSE)+
  scale_shape_manual(values=c(1,16)) + 
  scale_color_manual(values=c("grey", "black"))+
  geom_jitter(aes(color=group,shape=sig), width = 0.2, height = 0, show.legend = FALSE)+
  scale_x_discrete(labels=c("g1" = "Constance as P1/P2", "g2" = "C. acrinasus in P1 or P2"))+
  geom_signif(comparisons=list(c("g1", "g2")), annotations="***",
              y_position = 0.07, tip_length = 0.02, vjust=0.4)+
  theme_classic()
grid.arrange(dd,                              # box plot and scatter plot
             ncol = 1, nrow = 1)

#plots for FIGURE 2
tiff("Figures/Figure2_2.tiff", height=16, width=10, units="in", res=300, compression="lzw")
grid.arrange(grid.arrange(bb,aa,cc,dd, layout_matrix = rbind(c(1,1),c(4,NA), c(2,NA), c(3,NA))))
dev.off()

#test stats:
####stats: a)
wilcox.test(full_no047_nopel_wartmanii$Dstatistic, pelagic_nodup_wartmanii$Dstatistic, alternative = "less")
#W = 224, p-value < 2.2e-16

wilcox.test(full_no047_nopel_macrophthalmus$Dstatistic, pelagic_nodup_macrophthalmus$Dstatistic, alternative = "less")
#W = 3281, p-value = 0.0001391

wilcox.test(full_no047_nopel_arenicolus$Dstatistic, pelagic_nodup_arenicolus$Dstatistic, alternative = "less")
#W = 6173, p-value = 0.9172

####stats: b)
wilcox.test(full_con_nodups$Dstatistic, full_acr_nodups2$Dstatistic, alternative = "less")
#W = 11566, p-value < 2.2e-16

####stats: c)
wilcox.test(full_nopelagic_nonobilis$Dstatistic, prof_nodup_L$Dstatistic, alternative = "less")
#W = 70, p-value = 7.112e-12

####stats: d)
wilcox.test(full_noacrinasus_noprof$Dstatistic, prof_nodup_T$Dstatistic, alternative = "less")
#W = 247032, p-value = 0.9999

####numbers of topologies: a)
nrow(full_no047_nopel_wartmanii)
#144
nrow(pelagic_nodup_wartmanii)
#65
nrow(full_no047_nopel_macrophthalmus)
#99
nrow(pelagic_nodup_macrophthalmus)
#95
nrow(full_no047_nopel_arenicolus)
#102
nrow(pelagic_nodup_arenicolus)
#109

####numbers of topologies: b)
nrow(full_con_nodups)
#456
nrow(full_acr_nodups2)
#136

####numbers of topologies: c)
nrow(full_nopelagic_nonobilis)
#105
nrow(prof_nodup_L)
#21

####numbers of topologies: d)
nrow(full_noacrinasus_noprof)
#876
nrow(prof_nodup_T)
#502
