#set directory and load background file
setwd("/Users/rishide-kayne/Dropbox/RishiMAC/99_reanalysis/")
background <- read.csv("reviewer_responses/background_2021_99_resub.csv", header = T)

library(car)

background_filt <- background[-grep("missing", background$standard_length),]
background_filt2 <- background_filt[-grep("missing", background_filt$gill_raker_count),]

str(background_filt2)

background_filt2$symb <- c()
for (i in 1:nrow(background_filt2)){
  if(background_filt2$ecomorph_mod[i] == "balchen"){
    background_filt2$symb[i] <- 22
  }
  if(background_filt2$ecomorph_mod[i] == "albeli"){
    background_filt2$symb[i] <- 21
  }
  if(background_filt2$ecomorph_mod[i] == "felchen"){
    background_filt2$symb[i] <- 23
  }
  if(background_filt2$ecomorph_mod[i] == "large_pelagic"){
    background_filt2$symb[i] <- 8
  }
  if(background_filt2$ecomorph_mod[i] == "pelagic_profundal"){
    background_filt2$symb[i] <- 24
  }
  if(background_filt2$ecomorph_mod[i] == "benthic_profundal"){
    background_filt2$symb[i] <- 25
  }
}

background_filt2$standard_length <- as.numeric(background_filt2$standard_length)

background_filt2$gill_raker_count <- as.numeric(background_filt2$gill_raker_count)

par(mfrow=c(1,1))
plot(background_filt2$standard_length, background_filt2$gill_raker_count, 
     col = as.character(background_filt2$lake_col), pch = as.integer(background_filt2$symb), lwd = 2, ylab = "Gill raker count", xlab = "standard length (mm)", cex = 1.1)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Lucerne", "Biel", "Neuchâtel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1.1, cex = 0.75, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
test <- lm(background_filt2$standard_length ~ background_filt2$gill_raker_count)
summary(test)
abline(lm(background_filt2$gill_raker_count ~ background_filt2$standard_length), lwd = 2, lty = 3, col = "black")
#Multiple R-squared:  0.2767,	Adjusted R-squared:  0.2684 
#F-statistic: 33.28 on 1 and 87 DF,  p-value: 1.199e-07

l <- subset(background_filt2, background_filt2$Lake == "Lucerne")
b <- subset(background_filt2, background_filt2$Lake == "Brienz")
t <- subset(background_filt2, background_filt2$Lake == "Thun")
bi <- subset(background_filt2, background_filt2$Lake == "Biel")
n <- subset(background_filt2, background_filt2$Lake == "Neuenburg")
z <- subset(background_filt2, background_filt2$Lake == "Zurich")
w <- subset(background_filt2, background_filt2$Lake == "Walen")
c <- subset(background_filt2, background_filt2$Lake == "Constance")

par(mfrow=c(1,2))
boxplot(l$gill_raker_count, b$gill_raker_count, t$gill_raker_count, bi$gill_raker_count, n$gill_raker_count, z$gill_raker_count, w$gill_raker_count, c$gill_raker_count)
boxplot(l$standard_length, b$standard_length, t$standard_length, bi$standard_length, n$standard_length, z$standard_length, w$standard_length, c$standard_length)

par(mfrow=c(4,2))
par(mar=c(2,3,3,1))

test <- lm(l$standard_length ~ l$gill_raker_count)
plot(l$standard_length, l$gill_raker_count, col = l$lake_col, pch = l$symb)
summary(test)

test <- lm(c$standard_length ~ c$gill_raker_count)
plot(c$standard_length, c$gill_raker_count, col = c$lake_col, pch = c$symb)
summary(test)

test <- lm(b$standard_length ~ b$gill_raker_count)
plot(b$standard_length, b$gill_raker_count, col = b$lake_col, pch = b$symb)
summary(test)

test <- lm(t$standard_length ~ t$gill_raker_count)
plot(t$standard_length, t$gill_raker_count, col = t$lake_col, pch = t$symb)
summary(test)

test <- lm(bi$standard_length ~ bi$gill_raker_count)
plot(bi$standard_length, bi$gill_raker_count, col = bi$lake_col, pch = bi$symb)
summary(test)

test <- lm(n$standard_length ~ n$gill_raker_count)
plot(n$standard_length, n$gill_raker_count, col = n$lake_col, pch = n$symb)
summary(test)

test <- lm(z$standard_length ~ z$gill_raker_count)
plot(z$standard_length, z$gill_raker_count, col = z$lake_col, pch = z$symb)
summary(test)

test <- lm(w$standard_length ~ w$gill_raker_count)
plot(w$standard_length, w$gill_raker_count, col = w$lake_col, pch = w$symb)
summary(test)

l_bal <- subset(l, l$ecomorph_mod == "balchen")
l_alb <- subset(l, l$ecomorph_mod == "albeli")
l_fel <- subset(l, l$ecomorph_mod == "felchen")
l_lp <- subset(l, l$ecomorph_mod == "large_pelagic")
l_bp <- subset(l, l$ecomorph_mod == "benthic_prof")
l_pp <- subset(l, l$ecomorph_mod == "pelagic_profundal")

b_bal <- subset(b, b$ecomorph_mod == "balchen")
b_alb <- subset(b, b$ecomorph_mod == "albeli")
b_fel <- subset(b, b$ecomorph_mod == "felchen")
b_lp <- subset(b, b$ecomorph_mod == "large_pelagic")
b_bp <- subset(b, b$ecomorph_mod == "benthic_prof")
b_pp <- subset(b, b$ecomorph_mod == "pelagic_profundal")

t_bal <- subset(t, t$ecomorph_mod == "balchen")
t_alb <- subset(t, t$ecomorph_mod == "albeli")
t_fel <- subset(t, t$ecomorph_mod == "felchen")
t_lp <- subset(t, t$ecomorph_mod == "large_pelagic")
t_bp <- subset(t, t$ecomorph_mod == "benthic_prof")
t_pp <- subset(t, t$ecomorph_mod == "pelagic_profundal")

boxplot(l_bal$gill_raker_count,l_fel$gill_raker_count,l_alb$gill_raker_count,l_lp$gill_raker_count,l_bp$gill_raker_count,l_pp$gill_raker_count,
        b_bal$gill_raker_count,b_fel$gill_raker_count,b_alb$gill_raker_count,b_lp$gill_raker_count,b_bp$gill_raker_count,b_pp$gill_raker_count,
        t_bal$gill_raker_count,t_fel$gill_raker_count,t_alb$gill_raker_count,t_lp$gill_raker_count,t_bp$gill_raker_count,t_pp$gill_raker_count)

#with stat test of albeli/balchen overlap
par(mfrow=c(1,1))
plot(background_filt2$standard_length, background_filt2$gill_raker_count, 
     col = as.character(background_filt2$lake_col), pch = as.integer(background_filt2$symb), lwd = 2, ylab = "Gill raker count", xlab = "Standard length (mm)", cex = 1.1)
test <- lm(background_filt2$standard_length ~ background_filt2$gill_raker_count)
legend('topright', 
       legend = c("Brienz", "Thun", "Contance", "Walen", "Zurich", "Lucerne", "Biel", "Neuchâtel", "albeli", "felchen", 'balchen', "large pelagic", "benthic profundal", "pelagic profundal"), 
       col = c("chartreuse2", "chartreuse4", "goldenrod", "mediumpurple1", "purple4", "#ff4489", "paleturquoise2", "turquoise4", "black", "black", "black", "black", "black", "black"), 
       pt.cex = 1.1, cex = 0.55, pch = c(16, 16, 16, 16, 16, 16, 16, 16, 21, 23, 22, 8, 25, 24))
summary(test)
abline(lm(background_filt2$gill_raker_count ~ background_filt2$standard_length), lwd = 2, lty = 3, col = "black")

###ALL####
#make centroids
albeli <- subset(background_filt2, background_filt2$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(background_filt2, background_filt2$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)
#dataEllipse(albeli_x_vals, albeli_y_vals, levels=c(0.95), fill=TRUE, fill.alpha=0.1, plot.points=FALSE, col = "chartreuse2")
#dataEllipse(balchen_x_vals, balchen_y_vals, levels=c(0.95), fill=TRUE, fill.alpha=0.1, plot.points=FALSE, col = "chartreuse2")

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "grey", lty = 2)
#ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "grey", lty = 2, alpha = 0.75)

ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "grey", lty = 2)
#ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "grey", lty = 2, alpha = 0.75)

boxplot(albeli_x_vals, balchen_x_vals, ylab = "standard length")
wilcox.test(albeli_x_vals, balchen_x_vals)
#W = 33, p-value = 2.457e-07
#albeli n=19, balchen n=30

boxplot(albeli_y_vals, balchen_y_vals, ylab = "gill raker count")
wilcox.test(albeli_y_vals, balchen_y_vals)
#W = 533.5, p-value = 3.351e-07

hist(albeli_x_vals)
hist(balchen_x_vals)
hist(albeli_y_vals)
hist(balchen_y_vals)

###BRIENZ####
#make centroids
albeli <- subset(b, b$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(b, b$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)
#dataEllipse(albeli_x_vals, albeli_y_vals, levels=c(0.95), fill=TRUE, fill.alpha=0.1, plot.points=FALSE, col = "chartreuse2")
#dataEllipse(balchen_x_vals, balchen_y_vals, levels=c(0.95), fill=TRUE, fill.alpha=0.1, plot.points=FALSE, col = "chartreuse2")

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "chartreuse2", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "chartreuse2", lty = 2)

###WALEN####
albeli <- subset(w, w$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(w, w$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "mediumpurple1", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "mediumpurple1", lty = 2)

###NEUCHATEL####
albeli <- subset(n, n$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(n, n$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "turquoise4", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "turquoise4", lty = 2)

###LUCERENE####
albeli <- subset(l, l$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(l, l$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "#ff4489", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "#ff4489", lty = 2)

########################################################

###THUN####
albeli <- subset(t, t$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(t, t$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "chartreuse4", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "chartreuse4", lty = 2)

###BIEL####
albeli <- subset(bi, bi$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(bi, bi$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "paleturquoise2", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "paleturquoise2", lty = 2)

###SURICH####
albeli <- subset(z, z$ecomorph_mod == "albeli")
albeli_x <- (mean(albeli$standard_length))
albeli_y <- (mean(albeli$gill_raker_count))
points(albeli_x, albeli_y, pch = 16, cex = 2)

balchen <- subset(z, z$ecomorph_mod == "balchen")
balchen_x <- (mean(balchen$standard_length))
balchen_y <- (mean(balchen$gill_raker_count))
points(balchen_x, balchen_y, pch = 15, cex = 2)

#now plot ellipses
albeli_x_vals <- (albeli$standard_length)
albeli_y_vals <- (albeli$gill_raker_count)

albeli_df <- as.data.frame(cbind(albeli_x_vals, albeli_y_vals))
mu_albeli <- with(albeli_df, c(mean(albeli_x_vals), mean(albeli_y_vals)) )
sigma_albeli <- var(albeli_df)  # returns a variance-covariance matrix.

balchen_x_vals <- (balchen$standard_length)
balchen_y_vals <- (balchen$gill_raker_count)

balchen_df <- as.data.frame(cbind(balchen_x_vals, balchen_y_vals))
mu_balchen <- with(balchen_df, c(mean(balchen_x_vals), mean(balchen_y_vals)) )
sigma_balchen <- var(balchen_df)  # returns a variance-covariance matrix.


require(mixtools)

ellipse(mu_albeli, sigma_albeli, newplot = FALSE, col = "purple4", lty = 2)
ellipse(mu_balchen, sigma_balchen, newplot = FALSE, col = "purple4", lty = 2)


