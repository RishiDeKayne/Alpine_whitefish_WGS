###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	5.  GWAS
#   - 5.1 gill raker count
#		- 5.2 sex
#		- 5.3 standard length

#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)

################# 5. 1 GWAS gill raker count ############################### 
#new input has all outgroups and 06,070, and 071 removed since no gill raker count
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
emmax <- read.csv("GWAS/emmax_GR_SNP_pvalue_output_not0.05.txt", header = FALSE, sep = "\t")

#convert p value to -log10
emmax$log <- -log10(emmax$V3)

#get chromosome names and make vector
chroms <- read.csv(file = "Parallel/chrom_names.txt", header = FALSE)
chromlist <- as.list(as.character(chroms$V1))
b <- as.vector(chromlist)

a <- b[-c(22, 28, 32, 35, 38, 40)]

#would use pdf but too many points so causes crashes when viewing
jpeg(file="GWAS/test_GR.jpg",width=2000, height=600)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
#lengthsdf <- read.csv(file = "../Chap3/Fst/parallel_20_06_22/wtdbg2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- read.csv(file = "Parallel/wtdbg2ChromosomeLenghts.txt", header = FALSE)

lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
l <- as.vector(lengthsdf$V1)

#remove the absent chromosomes
l <- l[-c(22, 28, 32, 35, 38, 40)]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:34),nrow=1,ncol=34,byrow=TRUE),width=c(l))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("grey68", "grey50")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

#hard cuttoff is:
hard_cut <- -log10(0.05/9120498)

#LD cut is:
LD_cut <- -log10(0.05/4536915)

#extract snps with association over 7
outpoints_GR <- subset(emmax, emmax$log > LD_cut)

#and collect top % quantile of points
onepercent<-ecdf(emmax$log)
#then make new df with old data and the new p-value for each fst
onepercent_cutoff<-cbind(emmax,FST_p=(1-onepercent(emmax$log)))

onlyonepercent<-onepercent_cutoff[onepercent_cutoff$FST_p<0.0001,]

#now plot for each lake
#brienz
#for (i in 1:length(a)){
for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot

  plot(xFst,yFst,pch=".",cex=.5,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,10),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_GR$V1)){
    oneoutpoints_xFst <- outpoints_GR$V2[as.character(outpoints_GR$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_GR$log[as.character(outpoints_GR$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=1.5)}
  if (i==1){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
  rect(min(xFst), 0, max(xFst), -log10(0.05), col = "powderblue", border = NULL, lwd = 0.001)
}

dev.off()

emmax_GR <- emmax

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
tiff("GWAS/gill_raker_outlier_peak_freqs_by_ecomorph.tiff", height=5, width=7, units="in", res=300, compression="lzw")
par(mar=c(4,5,4,5)+0.1)
plot(full_snp_comparison$order_val, full_snp_comparison$V5, ylab = "reference allele frequency", xlab = "ecomorph", lwd = 2, xaxt='n', pch=c(25, 22, 23, 8, 21, 24), col = "black", ylim = c(0,1), cex = 1)
axis(1, at = c(1,2,3,4,5,6), labels=c("Benthic profundal", "Balchen", "Felchen", "Large pelagic", "Albeli", "Pelagic profundal"), las=1, cex = 1, cex.axis=0.7)
lines(full_snp_comparison$order_val,full_snp_comparison$V5, col = "grey50", lty = 3)

#and add the gill raker points
par(new = T)
plot(c(1,2,3,4,5,6), c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean), axes=F, xlab=NA, ylab=NA, cex=1, col = "orangered4", ylim = c(15, 40), pch=c(25, 22, 23, 8, 21, 24), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'mean gill raker count', col = "orangered4", cex = 1)
dev.off()

GRCs <- c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean)
AFs <- as.vector(full_snp_comparison$V5)
#now plot one against the other
tiff("GWAS/gill_raker_outlier_peak_freqs_by_ecomorph.tiff", height=5, width=7, units="in", res=300, compression="lzw")
par(mar=c(4,5,4,5)+0.1)
plot(GRCs, AFs, ylab = "reference allele frequency", xlab = "gill-raker count", lwd = 2, pch=c(25, 22, 23, 8, 21, 24), col = "black", ylim = c(0,1), cex = 1)
lines(GRCs,AFs, col = "grey50", lty = 3)
grc_lm <- as.data.frame(cbind(GRCs, AFs))

grc_lm$GRCs <- as.numeric(as.character(grc_lm$GRCs))
grc_lm$AFs <- as.numeric(as.character(grc_lm$AFs))

summary(lm(grc_lm$GRCs~grc_lm$AFs))
abline((lm(grc_lm$AFs~grc_lm$GRCs)), col = "red", lty = 2)

#and add the gill raker points
par(new = T)
plot(c(1,2,3,4,5,6), c(BP_gr_mean, balchen_gr_mean, felchen_gr_mean, LP_gr_mean, albeli_gr_mean, PP_gr_mean), axes=F, xlab=NA, ylab=NA, cex=1, col = "orangered4", ylim = c(15, 40), pch=c(25, 22, 23, 8, 21, 24), lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'mean gill raker count', col = "orangered4", cex = 1)
dev.off()

#now get % varience explained by most significant snp:
#PGA_scaffold22__199_contigs__length_52020451  30197713
#p-value = 7.962169e-09

#PGA_scaffold22__199_contigs__length_52020451	30197713	.	3.654445662	0.5726741251	7.962169472e-09

#N - sample size
#se_beta - standard error of effect size for the genetic variant of interest
#beta - effect size for the genetic variant of interest
#MAF - minor allele frequency for the genetic variant of interest

N <- 90
#then from the emmax output:
se_beta <- 0.5726813521
beta <- 3.654445661
  
#and from a vcftools command:
MAF <- 0.25

#then from : https://doi.org/10.1371/journal.pone.0120758.s001
#we get this equation for percentage variance explained:
PVE <- ((2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF)))
PVE


#and now lake-specific values:
luzern_freq <- read.csv("GWAS/luzern.frq.out", sep = "\t", head = FALSE)
luzern_freq$lake <- "L"
thun_freq <- read.csv("GWAS/thun.frq.out", sep = "\t", head = FALSE)
thun_freq$lake <- "T"
brienz_freq <-read.csv("GWAS/brienz.frq.out", sep = "\t", head = FALSE)
brienz_freq$lake <- "Br"
brienz_freq$V1 <- gsub("_B.frq", "_Br.frq", brienz_freq$V1)
zurich_freq <-read.csv("GWAS/zurich.frq.out", sep = "\t", head = FALSE)
zurich_freq$lake <- "Z"
walen_freq <- read.csv("GWAS/walen.frq.out", sep = "\t", head = FALSE)
walen_freq$lake <- "W"
biel_freq <- read.csv("GWAS/biel.frq.out", sep = "\t", head = FALSE)
biel_freq$lake <- "Bi"
neuchatel_freq <- read.csv("GWAS/neuchatel.frq.out", sep = "\t", head = FALSE)
neuchatel_freq$lake <- "N"
constance_freq <- read.csv("GWAS/constance.frq.out", sep = "\t", head = FALSE)
constance_freq$lake <- "C"

all.freqs <- rbind(luzern_freq, brienz_freq, thun_freq, zurich_freq, biel_freq, walen_freq, neuchatel_freq, constance_freq)
all.freqs$V6 <- gsub("C:", "", all.freqs$V6)

all.freqs$V7 <- gsub("A:", "", all.freqs$V7)

all.freqs$V1 <- gsub(".frq", "", all.freqs$V1)
all.freqs$V1 <- gsub("_prof", "", all.freqs$V1)
all.freqs$V1 <- gsub("_pel", "", all.freqs$V1)

library(tidyr)
new.all.freqs <- as.data.frame(separate(data = all.freqs, col = V1, into = c("ecomorph", "lake"), sep = "_"))
new.all.freqs$lake2 <- all.freqs$lake

new.all.freqs$meanGRC <- "empty"
new.all.freqs$col <- "empty"

lake_name_vec <- c()
e_name_vec <- c()
col_vec <- c()
for (i in 1:nrow(new.all.freqs)){
  if(new.all.freqs$lake2[i] == "L"){
    lake_name <- "Lucerne"
    col_plot <- "#ff4489"}
  if(new.all.freqs$lake2[i] == "T"){
    lake_name <- "Thun"
    col_plot <- "chartreuse4"}
  if(new.all.freqs$lake2[i] == "Br"){
    lake_name <- "Brienz"
    col_plot <- "chartreuse2"}
  if(new.all.freqs$lake2[i] == "Z"){
    lake_name <- "Zurich"
    col_plot <- "purple4"}
  if(new.all.freqs$lake2[i] == "W"){
    lake_name <- "Walen"
    col_plot <- "mediumpurple1"}
  if(new.all.freqs$lake2[i] == "Bi"){
    lake_name <- "Biel"
    col_plot <- "paleturquoise3"}
  if(new.all.freqs$lake2[i] == "N"){
    lake_name <- "Neuenburg"
    col_plot <- "turquoise4"}
  if(new.all.freqs$lake2[i] == "C"){
    lake_name <- "Constance"
    col_plot <- "goldenrod"}
  print(lake_name)
  lake_name_vec[i] <- lake_name
  print(col_plot)
  col_vec[i] <- col_plot
  if(new.all.freqs$ecomorph[i] == "albeli"){
    e_name <- "albeli"}
  if(new.all.freqs$ecomorph[i] == "balchen"){
    e_name <- "balchen"}
  if(new.all.freqs$ecomorph[i] == "felchen"){
    e_name <- "felchen"}
  if(new.all.freqs$ecomorph[i] == "large"){
    e_name <- "large_pelagic"}
  if(new.all.freqs$ecomorph[i] == "pel"){
    e_name <- "pelagic_profundal"}
  if(new.all.freqs$ecomorph[i] == "benth"){
    e_name <- "benthic_profundal"}
  print(e_name)
  e_name_vec[i] <- e_name
  
  datasub <- subset(background, as.character(background$Lake) == as.character(lake_name))
  datasub2 <- subset(datasub, as.character(datasub$ecomorph_mod) == as.character(e_name))
  datasub3 <- subset(datasub2, as.character(datasub2$gill_raker_count) != "missing")
  
  gill_vec <- as.vector(as.numeric(as.character(datasub3$gill_raker_count)))
  new.all.freqs$meanGRC[i] <- mean(gill_vec)
  new.all.freqs$col[i] <- col_plot
}

new.all.freqs$symb <- "empty"

for (i in 1:nrow(new.all.freqs)){
  if(new.all.freqs$ecomorph[i] == "albeli"){
    symb <- 21
  }
  if(new.all.freqs$ecomorph[i] == "balchen"){
    symb <- 22
  }
  if(new.all.freqs$ecomorph[i] == "felchen"){
    symb <- 23
  }
  if(new.all.freqs$ecomorph[i] == "large"){
    symb <- 8 
  }
  if(new.all.freqs$ecomorph[i] == "pel"){
    symb <- 24
  }
  if(new.all.freqs$ecomorph[i] == "benth"){
    symb <- 25
  }
  new.all.freqs$symb[i] <- as.character(symb)
}

plot(new.all.freqs$meanGRC, new.all.freqs$V6, pch = as.numeric(new.all.freqs$symb), col = as.character(new.all.freqs$col), lwd = 2, ylim = c(0,1))

luz <- subset(new.all.freqs, new.all.freqs$lake2 == "L")
luz <- luz[order(luz$meanGRC),]
lines(luz$meanGRC,luz$V6, col = "#ff4489", lty = 3)
thu <- subset(new.all.freqs, new.all.freqs$lake2 == "T")
thu <- thu[order(thu$meanGRC),]
lines(thu$meanGRC,thu$V6, col = "chartreuse4", lty = 3)
bri <- subset(new.all.freqs, new.all.freqs$lake2 == "Br")
bri <- bri[order(bri$meanGRC),]
lines(bri$meanGRC,bri$V6, col = "chartreuse2", lty = 3)
zur <- subset(new.all.freqs, new.all.freqs$lake2 == "Z")
zur <- zur[order(zur$meanGRC),]
lines(zur$meanGRC,zur$V6, col = "purple4", lty = 3)
wal <- subset(new.all.freqs, new.all.freqs$lake2 == "W")
wal <- wal[order(wal$meanGRC),]
lines(wal$meanGRC,wal$V6, col = "mediumpurple1", lty = 3)
bie <- subset(new.all.freqs, new.all.freqs$lake2 == "Bi")
bie <- bie[order(bie$meanGRC),]
lines(bie$meanGRC,bie$V6, col = "paleturquoise3", lty = 3)
neu <- subset(new.all.freqs, new.all.freqs$lake2 == "N")
neu <- neu[order(neu$meanGRC),]
lines(neu$meanGRC,neu$V6, col = "turquoise4", lty = 3)
con <- subset(new.all.freqs, new.all.freqs$lake2 == "C")
con <- con[order(con$meanGRC),]
lines(con$meanGRC,con$V6, col = "goldenrod", lty = 3)


new.all.freqs.noZalbeli <- new.all.freqs[-c(14), ]
plot(new.all.freqs.noZalbeli$meanGRC, new.all.freqs.noZalbeli$V6, 
     pch = as.numeric(new.all.freqs.noZalbeli$symb), 
     col = as.character(new.all.freqs.noZalbeli$col), 
     lwd = 2, 
     ylim = c(0,1),
     ylab = "Reference allele frequency",
     xlab = "Mean gill-raker count for ecomorph")

luz <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "L")
luz <- luz[order(luz$meanGRC),]
lines(luz$meanGRC,luz$V6, col = "#ff4489", lty = 3)
thu <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "T")
thu <- thu[order(thu$meanGRC),]
lines(thu$meanGRC,thu$V6, col = "chartreuse4", lty = 3)
bri <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "Br")
bri <- bri[order(bri$meanGRC),]
lines(bri$meanGRC,bri$V6, col = "chartreuse2", lty = 3)
zur <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "Z")
zur <- zur[order(zur$meanGRC),]
lines(zur$meanGRC,zur$V6, col = "purple4", lty = 3)
wal <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "W")
wal <- wal[order(wal$meanGRC),]
lines(wal$meanGRC,wal$V6, col = "mediumpurple1", lty = 3)
bie <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "Bi")
bie <- bie[order(bie$meanGRC),]
lines(bie$meanGRC,bie$V6, col = "paleturquoise3", lty = 3)
neu <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "N")
neu <- neu[order(neu$meanGRC),]
lines(neu$meanGRC,neu$V6, col = "turquoise4", lty = 3)
con <- subset(new.all.freqs.noZalbeli, new.all.freqs.noZalbeli$lake2 == "C")
con <- con[order(con$meanGRC),]
lines(con$meanGRC,con$V6, col = "goldenrod", lty = 3)


################# 5. 2 sex ############################### 
#new input has all outgroups and 06,070, and 071 removed since no gill raker count
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
emmax <- read.csv("GWAS/emmax_sex_SNP_pvalue_output_not0.05.txt", header = FALSE, sep = "\t")

#convert p value to -log10
emmax$log <- -log10(emmax$V3)

#get chromosome names and make vector
chroms <- read.csv(file = "Parallel/chrom_names.txt", header = FALSE)
chromlist <- as.list(as.character(chroms$V1))
b <- as.vector(chromlist)

a <- b[-c(22, 28, 32, 35, 38, 40)]

#would use pdf but too many points so causes crashes when viewing
jpeg(file="GWAS/test_sex.jpg",width=2000, height=600)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../Chap3/Fst/parallel_20_06_22/wtdbg2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
l <- as.vector(lengthsdf$V1)

#remove the absent chromosomes
l <- l[-c(22, 28, 32, 35, 38, 40)]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:34),nrow=1,ncol=34,byrow=TRUE),width=c(l))

#make colour vector to plot alternating chromosomes alternating colours
altcols <- c("grey68", "grey50")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

#extract snps with association over 7
outpoints_sex <- subset(emmax, emmax$log > 7.3)

#and collect top % quantile of points
onepercent<-ecdf(emmax$log)
#then make new df with old data and the new p-value for each fst
onepercent_cutoff<-cbind(emmax,FST_p=(1-onepercent(emmax$log)))

onlyonepercent<-onepercent_cutoff[onepercent_cutoff$FST_p<0.0001,]

#now plot for each lake
#brienz
#for (i in 1:length(a)){
for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,20),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_sex$V1)){
    oneoutpoints_xFst <- outpoints_sex$V2[as.character(outpoints_sex$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_sex$log[as.character(outpoints_sex$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=1.5)}
  if (i==1){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
  rect(min(xFst), 0, max(xFst), -log10(0.05), col = "powderblue", border = NULL, lwd = 0.001)}

dev.off()

N <- 90
#then from the emmax output:
se_beta <- 0.07238415936
beta <- 0.7410381949

#and from a vcftools command:
MAF <- 0.255556

#then from : https://doi.org/10.1371/journal.pone.0120758.s001
#we get this equation for percentage variance explained:
PVE_sex <- ((2*(beta^2)*MAF*(1-MAF))/(2*(beta^2)*MAF*(1-MAF)+(se_beta^2)*2*N*MAF*(1-MAF)))
PVE_sex

#############now both together

jpeg(file="GWAS/GR_and_sex_GWAS_largepoints.jpg",width=2000, height=1300)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "../Chap3/Fst/parallel_20_06_22/wtdbg2ChromosomeLenghts.txt", header = FALSE)
lengthsdf <- subset(lengthsdf, lengthsdf$V1 > 50)
l <- as.vector(lengthsdf$V1)

#remove the absent chromosomes
l <- l[-c(22, 28, 32, 35, 38, 40)]

#set layout parameters specifying 5 rows, one for each comarison and the widths of the chromosomes in each case
layout(matrix(seq(1:68),nrow=2,ncol=34,byrow=TRUE),width=c(l))

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax_GR$V2[as.character(emmax_GR$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax_GR$log[as.character(emmax_GR$V1)==as.character(a[i])]
  #then plot
  
  plot(xFst,yFst,pch=16,cex=1,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,10),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_GR$V1)){
    oneoutpoints_xFst <- outpoints_GR$V2[as.character(outpoints_GR$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_GR$log[as.character(outpoints_GR$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=2.5)}
  if (i==1){
    axis(side=2)
  }
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
  rect(min(xFst), 0, max(xFst), -log10(0.05), col = "powderblue", border = NULL, lwd = 0.001)
}

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot
  plot(xFst,yFst,pch=16,cex=1,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,20),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_sex$V1)){
    oneoutpoints_xFst <- outpoints_sex$V2[as.character(outpoints_sex$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_sex$log[as.character(outpoints_sex$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=2.5)}
  if (i==1){
    axis(side=2)
  }
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
  rect(min(xFst), 0, max(xFst), -log10(0.05), col = "powderblue", border = NULL, lwd = 0.001)}
dev.off()

##############################
# now run on cluster to plot all points 
setwd("/cluster/scratch/rdekayne/99_reanalysis/GWAS_poly/figs/")
#background <- read.csv("background_2021_99.csv", header = T)

################# 5. 1 GWAS gill raker count ############################### 
#new input has all outgroups and 06,070, and 071 removed since no gill raker count
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
emmax <- read.csv("../gill_raker/emmax_GR_SNP_pvalue_output_not1.txt", header = FALSE, sep = "\t")

#convert p value to -log10
emmax$log <- -log10(emmax$V3)

#get chromosome names and make vector
chroms <- read.csv(file = "chrom_names.txt", header = FALSE)
chromlist <- as.list(as.character(chroms$V1))
b <- as.vector(chromlist)

a <- b[-c(22, 28, 32, 35, 38, 40)]

#would use pdf but too many points so causes crashes when viewing
pdf(file="./GR_gwas.pdf",width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
altcols <- c("grey80", "grey60")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

#hard cuttoff is:
hard_cut <- -log10(0.05/9120498)

#LD cut is:
LD_cut <- -log10(0.05/4536915)

#extract snps with association over 7
outpoints_GR <- subset(emmax, emmax$log > LD_cut)

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot
  
  plot(xFst,yFst,pch=".",cex=.5,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,10),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_GR$V1)){
    oneoutpoints_xFst <- outpoints_GR$V2[as.character(outpoints_GR$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_GR$log[as.character(outpoints_GR$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=1.5)}
  if (i==1){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
  #rect(min(xFst), 0, max(xFst), -log10(0.05), col = "powderblue", border = NULL, lwd = 0.001)
}

dev.off()

emmax_GR <- emmax

################# 5. 2 sex ############################### 
#new input has all outgroups and 06,070, and 071 removed since no gill raker count
#THE FOLLOWING FILE IS TOO LARGE FOR GITHUB AND WILL BE ARCHIVED ELSEWHERE - PLEASE SEE THE PAPER FOR DETAILS
emmax <- read.csv("../sex/emmax_sex_SNP_pvalue_output_not1.txt", header = FALSE, sep = "\t")

#convert p value to -log10
emmax$log <- -log10(emmax$V3)

#get chromosome names and make vector
chroms <- read.csv(file = "./chrom_names.txt", header = FALSE)
chromlist <- as.list(as.character(chroms$V1))
b <- as.vector(chromlist)

a <- b[-c(22, 28, 32, 35, 38, 40)]

#would use pdf but too many points so causes crashes when viewing
pdf(file="./sex_gwas.pdf",width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "./wtdbg2ChromosomeLenghts.txt", header = FALSE)
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
altcols <- c("grey80", "grey60")
altcols <- rep(altcols, as.integer((length(a)+1)/2))

#extract snps with association over LD_cut
outpoints_sex <- subset(emmax, emmax$log > LD_cut)

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot
  plot(xFst,yFst,pch=".",cex=.5,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,20),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_sex$V1)){
    oneoutpoints_xFst <- outpoints_sex$V2[as.character(outpoints_sex$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_sex$log[as.character(outpoints_sex$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=1.5)}
  if (i==1){
    axis(side=2)
  }
  axis(side=1,at=c(0,max(xFst)),outer=F,labels=F,tick=T,lwd.ticks=0,lwd=5,col='snow4')
  #abline(h = -log10(0.01), col = "red")
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
}

dev.off()

#############now both together

pdf(file="./GR_and_sex_gwas.pdf",width=15, height=4)

par(mfrow=c(1,1),mar=c(0.01,0.01,0.01,0.01),oma=c(4,3,3,1))

#use the length of each chromosome to plot relative widths on the skyline plots
lengthsdf <- read.csv(file = "./wtdbg2ChromosomeLenghts.txt", header = FALSE)
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

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax_GR$V2[as.character(emmax_GR$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax_GR$log[as.character(emmax_GR$V1)==as.character(a[i])]
  #then plot
  
  plot(xFst,yFst,pch=16,cex=1,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,10),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_GR$V1)){
    oneoutpoints_xFst <- outpoints_GR$V2[as.character(outpoints_GR$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_GR$log[as.character(outpoints_GR$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=2.5)}
  if (i==1){
    axis(side=2)
  }
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
}

for (i in 1:34){
  #extract the x/bp position for each chrom
  xFst <- emmax$V2[as.character(emmax$V1)==as.character(a[i])]
  #extract the fst values
  yFst <- emmax$log[as.character(emmax$V1)==as.character(a[i])]
  #then plot
  plot(xFst,yFst,pch=16,cex=1,ylab="-log10(p-value)",xaxt="n",axes=F,ylim=c(0,20),col=altcols[i])
  if (as.character(a[i]) %in% as.character(outpoints_sex$V1)){
    oneoutpoints_xFst <- outpoints_sex$V2[as.character(outpoints_sex$V1)==as.character(a[i])]
    oneoutpoints_yFst <- outpoints_sex$log[as.character(outpoints_sex$V1)==as.character(a[i])]
    points(oneoutpoints_xFst,oneoutpoints_yFst,col="orangered4",pch=20,cex=2.5)}
  if (i==1){
    axis(side=2)
  }
  abline(h = hard_cut, col = "orangered4", lty = 1)
  abline(h = LD_cut, col = "orangered4", lty = 2)
}

dev.off()