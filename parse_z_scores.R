#test - convert z-scores

zscores <- read.csv("Dropbox/RishiMAC/99_reanalysis/reviewer_responses/dsuite/z_scores.txt", sep = "\t")

header <- zscores[,1:2]
             
scores <- zscores[,3:ncol(zscores)]

newscores <- as.data.frame(c())
na_count <- 0
for(i in 1:nrow(scores)){
  for(j in 1:ncol(scores)){
    if(scores[i,j] != 'NaN'){
      if(as.numeric(scores[i,j])<3){
        newscores[i,j] <- 0
      }
      if(as.numeric(scores[i,j])>=3){
        newscores[i,j] <- 0.99
      }
    }
    if(scores[i,j] == 'NaN'){
      newscores[i,j] <- "nan"
      na_count <- na_count+1
    }
  }
}

all_new <- cbind(header, newscores)

colnames(all_new) <- colnames(zscores)

write.table(all_new, "Dropbox/RishiMAC/99_reanalysis/reviewer_responses/dsuite/all_new_zscore.txt", quote=FALSE, row.names=FALSE, sep = "\t")

#significance cutoff
cell_count <- nrow(newscores)*ncol(newscores)
cell_count_nona <- cell_count-na_count
pvalue <- 0.01

corrected_pvalue <- pvalue/cell_count_nona
corrected_pvalue
Zscore <- qnorm(corrected_pvalue, lower.tail=F)

newscores_Zscore <- as.data.frame(c())
na_count <- 0
for(i in 1:nrow(scores)){
  for(j in 1:ncol(scores)){
    if(scores[i,j] != 'NaN'){
      if(as.numeric(scores[i,j])<=Zscore){
        newscores_Zscore[i,j] <- 0
      }
      if(as.numeric(scores[i,j])>Zscore){
        newscores_Zscore[i,j] <- 0.99
      }
    }
    if(scores[i,j] == 'NaN'){
      newscores_Zscore[i,j] <- "nan"
      na_count <- na_count+1
    }
  }
}

all_new_zscore_sig <- cbind(header, newscores_Zscore)

colnames(all_new_zscore_sig) <- colnames(zscores)

write.table(all_new_zscore_sig, "Dropbox/RishiMAC/99_reanalysis/reviewer_responses/dsuite/all_new_zscore_sig.txt", quote=FALSE, row.names=FALSE, sep = "\t")


#then get rid of header 'X' added by R: sed -i 's/X//g' all_new_zscore.txt 

