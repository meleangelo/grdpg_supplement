rm(list=ls())
library(xtable)

setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc/")

table_est <- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_est) <- c("n", "betahat[1]", "betahat[2]", "betahat_w[1]", "betahat_w[2]", "time")
# n=2000
load("mc_estimates_nsim1000_n2000_all.RData")
table_est[1,1] <- 2000
table_est[1,2:6] <- colMeans(estimates_grdpg)
# n=5000
load("mc_estimates_nsim1000_n5000_all.RData")
table_est[2,1] <- 5000
table_est[2,2:6] <- colMeans(estimates_grdpg)
# n=10000
load("mc_estimates_nsim1000_n10000_all.RData")
table_est[3,1] <- 10000
table_est[3,2:6] <- colMeans(estimates_grdpg)

# latex table for paper
xtable(table_est, digits = c(0,0,4,4,4,4,2))


table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,2]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


#### table bias for weighted estimates
table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,4]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


##################################################
##### monte carlo with  correlated covariates
##################################################
rm(list=ls())
library(xtable)

setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc/")

table_est <- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_est) <- c("n", "betahat[1]", "betahat[2]", "betahat_w[1]", "betahat_w[2]", "time")
# n=2000
load("mc_correl_estimates_nsim1000_n2000_all.RData")
table_est[1,1] <- 2000
table_est[1,2:6] <- colMeans(estimates_grdpg)
# n=5000
load("mc_correl_estimates_nsim1000_n5000_all.RData")
table_est[2,1] <- 5000
table_est[2,2:6] <- colMeans(estimates_grdpg)
# n=10000
load("mc_correl_estimates_nsim1000_n10000_all.RData")
table_est[3,1] <- 10000
table_est[3,2:6] <- colMeans(estimates_grdpg)

# latex table for paper
xtable(table_est, digits = c(0,0,4,4,4,4,2))


table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_correl_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_correl_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_correl_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,2]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


#### table bias for weighted estimates
table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_correl_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,3] - 0.5), na.rm = T)
table_bias[1,3] <- sd(estimates_grdpg[,3], na.rm = T) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,4] - 0.75), na.rm = T)
table_bias[1,5] <- sd(estimates_grdpg[,4], na.rm = T)  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5], na.rm = T)

# n=5000
load("mc_correl_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_correl_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,4]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))




####################################################
# monte carlo with unbalanced blocks and indep covariates
####################################################
rm(list=ls())
library(xtable)

setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc/")


table_est <- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_est) <- c("n", "betahat[1]", "betahat[2]", "betahat_w[1]", "betahat_w[2]", "time")
# n=2000
load("mc_unbalanced_estimates_nsim1000_n2000_all.RData")
table_est[1,1] <- 2000
table_est[1,2:6] <- colMeans(estimates_grdpg)
# n=5000
load("mc_unbalanced_estimates_nsim1000_n5000_all.RData")
table_est[2,1] <- 5000
table_est[2,2:6] <- colMeans(estimates_grdpg)
# n=10000
load("mc_unbalanced_estimates_nsim1000_n10000_all.RData")
table_est[3,1] <- 10000
table_est[3,2:6] <- colMeans(estimates_grdpg)

# latex table for paper
xtable(table_est, digits = c(0,0,4,4,4,4,2))


table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_unbalanced_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_unbalanced_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_unbalanced_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,2]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


#### table bias for weighted estimates
table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_unbalanced_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_unbalanced_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_unbalanced_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,4]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))





####################################################
# monte carlo with unbalanced blocks and unbalanced covariates
####################################################

rm(list=ls())
library(xtable)

setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc/")

table_est <- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_est) <- c("n", "betahat[1]", "betahat[2]", "betahat_w[1]", "betahat_w[2]", "time")
# n=2000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_est[1,1] <- 2000
table_est[1,2:6] <- colMeans(estimates_grdpg)
# n=5000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_est[2,1] <- 5000
table_est[2,2:6] <- colMeans(estimates_grdpg)
# n=10000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_est[3,1] <- 10000
table_est[3,2:6] <- colMeans(estimates_grdpg)

# latex table for paper
xtable(table_est, digits = c(0,0,4,4,4,4,2))


table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,1] - 0.5), na.rm = T)
table_bias[1,3] <- sd(estimates_grdpg[,1], na.rm = T) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,2] - 0.75), na.rm = T)
table_bias[1,5] <- sd(estimates_grdpg[,2], na.rm = T)  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5], na.rm = T)

# n=5000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,2]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


#### table bias for weighted estimates
table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,3] - 0.5), na.rm = T)
table_bias[1,3] <- sd(estimates_grdpg[,3], na.rm = T) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,4] - 0.75), na.rm = T)
table_bias[1,5] <- sd(estimates_grdpg[,4], na.rm = T)  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5], na.rm = T)

# n=5000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,4])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,3] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,3]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,4] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,4]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))





####################################################
# monte carlo with unbalanced blocks and correlated unbalanced covariates
####################################################

rm(list=ls())
library(xtable)

setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc/")

table_est <- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_est) <- c("n", "betahat[1]", "betahat[2]", "betahat_w[1]", "betahat_w[2]", "time")
# n=2000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_est[1,1] <- 2000
table_est[1,2:6] <- colMeans(estimates_grdpg)
# n=5000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_est[2,1] <- 5000
table_est[2,2:6] <- colMeans(estimates_grdpg)
# n=10000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_est[3,1] <- 10000
table_est[3,2:6] <- colMeans(estimates_grdpg)

# latex table for paper
xtable(table_est, digits = c(0,0,4,4,4,4,2))


table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[1,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[1,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5])

# n=5000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[2,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[2,5] <- sd(estimates_grdpg[,2])  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5])

# n=10000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,1] - 0.5))
table_bias[3,3] <- sd(estimates_grdpg[,1]) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,2] - 0.75))
table_bias[3,5] <- sd(estimates_grdpg[,2]) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5])

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


#### table bias for weighted estimates
table_bias<- data.frame(matrix(NA, ncol = 6, nrow = 3))
names(table_bias) <- c("n", "betahat[1]", "mcse[1]", "betahat[2]", "mcse[2]", "time")

# n=2000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n2000_all.RData")
table_bias[1,1] <- 2000
table_bias[1,2] <- mean(abs(estimates_grdpg[,3] - 0.5), na.rm = T)
table_bias[1,3] <- sd(estimates_grdpg[,3], na.rm = T) /sqrt(n)
table_bias[1,4] <- mean(abs(estimates_grdpg[,4] - 0.75), na.rm = T)
table_bias[1,5] <- sd(estimates_grdpg[,4], na.rm = T)  /sqrt(n)
table_bias[1,6] <- mean(estimates_grdpg[,5], na.rm = T)

# n=5000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n5000_all.RData")
table_bias[2,1] <- 5000
table_bias[2,2] <- mean(abs(estimates_grdpg[,3] - 0.5), na.rm = T)
table_bias[2,3] <- sd(estimates_grdpg[,3], na.rm = T) /sqrt(n)
table_bias[2,4] <- mean(abs(estimates_grdpg[,4] - 0.75), na.rm = T)
table_bias[2,5] <- sd(estimates_grdpg[,4], na.rm = T)  /sqrt(n)
table_bias[2,6] <- mean(estimates_grdpg[,5], na.rm = T)

# n=10000
load("mc_correl_unbalanced_unbalcov_estimates_nsim1000_n10000_all.RData")
table_bias[3,1] <- 10000
table_bias[3,2] <- mean(abs(estimates_grdpg[,3] - 0.5), na.rm = T)
table_bias[3,3] <- sd(estimates_grdpg[,3], na.rm = T) /sqrt(n)
table_bias[3,4] <- mean(abs(estimates_grdpg[,4] - 0.75), na.rm = T)
table_bias[3,5] <- sd(estimates_grdpg[,4], na.rm = T) /sqrt(n)
table_bias[3,6] <- mean(estimates_grdpg[,5], na.rm = T)

# latex table for paper
xtable(table_bias, digits = c(0,0,4,7,4,7,2))


