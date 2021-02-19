rm(list=ls())
library(xtable)
datafolder <- "/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/" # macbook 15
setwd(datafolder)

# initialize empty table for estimates
table_est <- data.frame(matrix(NA, nrow = 9, ncol = 4))

# first column is gender only
load("Harvard_elbow1_allyears_gender.RData")
table_est[1,1] <- betahat2
table_est[2,1] <- sd2
table_est[7,1] <- dim(A)[1]
table_est[8,1] <- model$G
table_est[9,1] <- dhat

# second column is role only
load("Harvard_elbow1_allyears_role.RData")
table_est[3,2] <- betahat2
table_est[4,2] <- sd2
table_est[7,2] <- dim(A)[1]
table_est[8,2] <- model$G
table_est[9,2] <- dhat


# third column is dorm only
load("Harvard_elbow1_allyears_dorm.RData")
table_est[5,3] <- betahat2
table_est[6,3] <- sd2
table_est[7,3] <- dim(A)[1]
table_est[8,3] <- model$G
table_est[9,3] <- dhat

# fourth column is all variable
load("Harvard_elbow1_allyears_all.RData")
table_est[1,4] <- betahat2[1]
table_est[2,4] <- sd2[1]
table_est[3,4] <- betahat2[2]
table_est[4,4] <- sd2[2]
table_est[5,4] <- betahat2[3]
table_est[6,4] <- sd2[3]
table_est[7,4] <- dim(A)[1]
table_est[8,4] <- model$G
table_est[9,4] <- dhat

#rownames(table_est) <- c("female", "", "student", "", "off-campus", "", "n", "K", "d" )

xtable(table_est, digits = 4)



