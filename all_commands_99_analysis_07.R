###This script contains commands run for the whole-genome analysis paper for whitefish De-Kayne et al. 

#for all bash analyses please look at all_commands_99_analysis.txt

#then the following analyses

#	7.  f4 statistics
#Pop1 – brienz
#Pop2 – neuchatel
#Pop3 – Luzern
#Pop4 – walen 

#set directory and load background file
setwd("/Users/rishidek/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("background_2021_99.csv", header = T)


#read in data
data <- read.csv("f4/f4_output.txt", sep = ";", header = TRUE)
data$number <- 12:1

#arrange in plot so topologies are split up
ylabel_vector_norm <- as.character(data$topology)
ylabel_vector <- rev(ylabel_vector_norm)
par(mar=c(5,12,4,1)+.1)
plot(data$f4, data$number, pch = 16, ylab = '', yaxt='n', xlab = "f4")
axis(2, at = c(1:12), labels=ylabel_vector, las=1, cex = 1)
abline(h=6.5, col = "red", lty = 3)

#now want to illustrate comparison of tree topologies so want same 4 taxon tree rearranged for each combo
lakepattern <- data[1:6,]
lakepattern$newplotnumber <- c(6:1)
lakepattern$colour <- "darkslategray4"
ecotypepattern <- data[7:12,]
ecotypepattern$newplotnumber <- c(6:1)
ecotypepattern$colour <- "goldenrod3"

laketop <- as.character(lakepattern$topology)
lakef4 <- as.character(lakepattern$f4)
ecotop <- as.character(ecotypepattern$topology)
ecof4 <- as.character(ecotypepattern$f4)

sidebyside <- as.data.frame(cbind(laketop, lakef4, ecotop, ecof4))
str(sidebyside)

sidebyside$lakef4 <- as.numeric(as.character(sidebyside$lakef4))
sidebyside$ecof4 <- as.numeric(as.character(sidebyside$ecof4))

#calculate the difference between topologies
sidebyside$difference_in_value <- ((as.numeric(sidebyside$ecof4))-(as.numeric(sidebyside$lakef4)))

#now actually make plot
labels_norm <- labels_norm <- c("Brienz-Neuchâtel", "Brienz-Lucerne", "Brienz-Walen", "Neuchâtel-Lucerne", "Neuchâtel-Walen", "Lucerne-Walen")
labels_rev <- rev(labels_norm)

newplotdf <- rbind(lakepattern, ecotypepattern)
par(mar=c(5,8,4,1)+.1)

plot(newplotdf$f4, newplotdf$newplotnumber, pch = 16, ylab = '', yaxt='n', xlab = "f4", col = newplotdf$colour, ylim = c(0.5, 6.5))
axis(2, at = c(1:6), labels=labels_rev, las=1, cex = 1)
abline(h=0.5, col = "black", lty = 3)
abline(h=1.5, col = "black", lty = 3)
abline(h=2.5, col = "black", lty = 3)
abline(h=3.5, col = "black", lty = 3)
abline(h=4.5, col = "black", lty = 3)
abline(h=5.5, col = "black", lty = 3)
abline(h=6.5, col = "black", lty = 3)
mtext("(Big,Small);(Big,Small)", line = 1, side = 3, adj = 0.15)
mtext("(Big,Big);(Small,Small)", line = 1, side = 3, adj = 0.85)

newplotdf$newcol <- c()
for (i in 1:nrow(newplotdf)){
  if(newplotdf$colour[i] == "darkslategray4"){
    newplotdf$newcol[i] <- "grey"
  }
  else{
    newplotdf$newcol[i] <- "black"
  }
}
plot(newplotdf$f4, newplotdf$newplotnumber, pch = 16, ylab = '', yaxt='n', xlab = "f4", col = newplotdf$newcol, ylim = c(0.5, 6.5))
axis(2, at = c(1:6), labels=labels_rev, las=1, cex = 1)
abline(h=0.5, col = "black", lty = 3)
abline(h=1.5, col = "black", lty = 3)
abline(h=2.5, col = "black", lty = 3)
abline(h=3.5, col = "black", lty = 3)
abline(h=4.5, col = "black", lty = 3)
abline(h=5.5, col = "black", lty = 3)
abline(h=6.5, col = "black", lty = 3)
mtext("(Balchen,Albeli);(Balchen,Albeli)", line = 1, side = 3, adj = 0.15)
mtext("(Balchen,Balchen);(Albeli,Albeli)", line = 1, side = 3, adj = 0.85)


#raw data
#Pop1 – brienz
#Pop2 – neuchatel
#Pop3 – Luzern
#Pop4 – walen 

#(LGR1,LGR3),(HGR1,HGR3);0.02380
#(LGR1,HGR1),(LGR3,HGR3);0.01021
#(LGR1,LGR4),(HGR1,HGR4);0.02480
#(LGR1,HGR1),(LGR4,HGR4);0.00868
#(LGR2,LGR3),(HGR2,HGR3);0.02067
#(LGR2,HGR2),(LGR3,HGR3);0.00732
#(LGR2,LGR4),(HGR2,HGR4);0.01697
#(LGR2,HGR2),(LGR4,HGR4);0.00560
#(LGR1,LGR2),(HGR1,HGR2);0.02063
#(LGR3,LGR4),(HGR3,HGR4);0.02035
#(LGR3,HGR3),(LGR4,HGR4);0.01272
#(LGR1,HGR1),(LGR2,HGR2);0.00325

