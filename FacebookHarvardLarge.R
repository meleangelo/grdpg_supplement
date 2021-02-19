# Facebook analysis, Harvard University

rm(list = ls()) # clean memory

# required packages
require(dplyr)
require(mclust)
require(grdpg)
require(R.matlab)
require(igraph)

#### Data Preparation
# change next line to the folder containing the data
#datafolder <- "D:/grdpg_simulations/facebook100/" # research pc
datafolder <- "/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/" # macbook 15
# Harvard University data
file <- paste(datafolder, "Harvard1.mat", sep = "")


## Load data 
data.network <- readMat(file)

## Adjacency matrix
A <- as.matrix(data.network$A)

## Nodes information
nodeinfo <- data.frame(data.network$local.info)
vblnames <- c("role", "gender", "major", "minor", "dorm", "year",  "highschool")
names(nodeinfo) <- vblnames

# count size
original_data_size <-dim(A)[1]

## descriptive statistics
g <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
# avg degree
degree_orig <- igraph::degree(g)
avgdegree_orig <- mean(degree_orig)
# clustering
clust_orig <- igraph::transitivity(g, type = "undirected")
# % female, % students, % off-campus
female_pct_orig <- sum(nodeinfo$gender==2)/original_data_size
student_pct_orig <- sum(nodeinfo$role==1)/original_data_size
offcampus_pct_orig <- sum(nodeinfo$dorm==0)/original_data_size


## Drop if gender is missing
missing.gender <- which(nodeinfo$gender==0)
A <- A[-missing.gender,-missing.gender]
nodeinfo <- nodeinfo[-missing.gender,]

# count dimension
network_gender_size <- dim(A)[1]

# avg degree
g <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
degree_gender <- igraph::degree(g)
avgdegree_gender <- mean(degree_gender)
# clustering
clust_gender <- igraph::transitivity(g, type = "undirected")
# % female, % students, % off-campus
female_pct_gender <- sum(nodeinfo$gender==2)/network_gender_size
student_pct_gender <- sum(nodeinfo$role==1)/network_gender_size
offcampus_pct_gender <- sum(nodeinfo$dorm==0)/network_gender_size



# keep only largest connected component
components <- igraph::clusters(g)
biggest_cluster_id <- which.max(components$csize)
vert_ids <- igraph::V(g)[components$membership == biggest_cluster_id]
vert_largecomp <- as.numeric(vert_ids)
# take adjacency matrix and node info for largest component
A <- A[vert_largecomp,vert_largecomp]
nodeinfo <- nodeinfo[vert_largecomp,]

# count dimension
lcc_size <- dim(A)[1]


# avg degree
g <- igraph::graph_from_adjacency_matrix(A, mode="undirected")
degree_lcc <- igraph::degree(g)
avgdegree_lcc <- mean(degree_lcc)
# clustering
clust_lcc <- igraph::transitivity(g, type = "undirected")
# % female, % students, % off-campus
female_pct_lcc <- sum(nodeinfo$gender==2)/lcc_size
student_pct_lcc <- sum(nodeinfo$role==1)/lcc_size
offcampus_pct_lcc <- sum(nodeinfo$dorm==0)/lcc_size





## TABLES WITH DESCRIPTIVE STATS
# original data
original_data_size
avgdegree_orig
clust_orig
female_pct_orig
student_pct_orig
offcampus_pct_orig
# final network
network_gender_size
avgdegree_gender
clust_gender
female_pct_gender
student_pct_gender
offcampus_pct_gender
# largest connetd component
lcc_size
avgdegree_lcc
clust_lcc
female_pct_lcc
student_pct_lcc
offcampus_pct_lcc



## REGULARIZATION
# compute average degree
avgdegree <- sum(A)/(dim(A)[1])
A_orig <- A                               # keep original adjacency matrix
diag(A) <- diag(A)+ rowSums(A)/dim(A)[1]  # diagonal is augmented with avg degree/n
A <- A + avgdegree/dim(A)[1]              # Levina et al regularization for sparse graphs



## ANALYSIS WITH SINGLE COVARIATE: gender

cov <- 2
# modify the role variable to have only 2 values (student or not student)
# otherwise there are too few observations in some blocks
nodeinfo[nodeinfo$role>1,1] <- 2
# rewrite dorm as in campus/off-campus
nodeinfo[nodeinfo$dorm>0,5] <- 2
nodeinfo[nodeinfo$dorm==0,5] <- 1


covariates <- matrix(nodeinfo$gender, ncol = 1)


#### GRDPG with Covariates for Facebook Data
dmax <- 40 
G <- seq(from = 2, to = 50, by = prod(cov))

## ASE
#set.seed(1977)
set.seed(197)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension, use first elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[1]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile)
print("gender , all years, Harvard, 1st elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow1_allyears_gender.Rdata", sep = "")
save.image(file = file_estimates )


### choose d form second elbow
## Choose embed dimension, use second elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[2]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("gender ,  all years, Harvard, 2nd elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow2_allyears_gender.Rdata", sep = "")
save.image(file = file_estimates )





## ANALYSIS WITH SINGLE COVARIATE: role

cov <- 2
covariates <- matrix(nodeinfo$role, ncol = 1)


#### GRDPG with Covariates for Facebook Data
dmax <- 40 
G <- seq(from = 2, to = 50, by = prod(cov))

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension, use first elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[1]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("role , all years, Harvard, 1st elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow1_allyears_role.Rdata", sep = "")
save.image(file = file_estimates )


### choose d form second elbow
## Choose embed dimension, use second elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[2]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("role , all years, Harvard, 2nd elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow2_allyears_role.Rdata", sep = "")
save.image(file = file_estimates )



## ANALYSIS WITH SINGLE COVARIATE: off-campus

cov <- 2

covariates <- matrix(nodeinfo$dorm, ncol = 1)


#### GRDPG with Covariates for Facebook Data
dmax <- 40 
G <- seq(from = 2, to = 50, by = prod(cov))

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension, use first elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[1]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("off-campus , all years, Harvard, 1st elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow1_allyears_dorm.Rdata", sep = "")
save.image(file = file_estimates )


### choose d form second elbow
## Choose embed dimension, use second elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[2]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("offcampus , all years, Harvard, 2nd elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow2_allyears_dorm.Rdata", sep = "")
save.image(file = file_estimates )




## ANALYSIS WITH MULTIPLE COVARIATES (3 binary covariates)


cov <-c(2,2,2) 

covariates <- matrix(as.matrix(nodeinfo[,c(2,1,5)]), ncol = 3)


#### GRDPG with Covariates for Facebook Data
dmax <- 40 
G <- seq(from = prod(cov), to = prod(cov)*10, by = prod(cov))

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension, use first elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[1]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("gender + role + dorm (binary), all years , Harvard fixed blocks, 1st elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow1_allyears_all.Rdata", sep = "")
save.image(file = file_estimates )


### choose d form second elbow
## Choose embed dimension, use second elbow (Zhu and Ghodsi 2006)
dimselect(s)
dhat <- dimselect(s)$elbow[2]+1

## Estimate latent position
Yhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Yhat, G)
muhats <- model$parameters$mean
BXhat <- logit(BFcheck(t(muhats) %*% Ipq %*% muhats))

## Estimate beta
clusters_cov <- getClusters(data.frame(model$z))
covariates_block <- getBlockCovariates(covariates, clusters_cov)
check = "BF"


# Use weighted estimator
result2 <- estimatebeta_WA(Yhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
varweights <- Map('*', result2$pis, result2$pis)
var2 <- sapply(Map('*',result2$sd2s,varweights), sum)
sd2 <- sqrt(var2)/dim(A)[1]
betahat2 
betahat2_unbiased
sd2
sinkfile <- paste(datafolder, "table.txt", sep = "")
sink(sinkfile, append = TRUE)
print("gender + role + dorm (binary), all years, Harvard, 2nd elbow")
betahat2 
betahat2_unbiased  # naive bias correction
sd2                # naive standard error
sink()

file_estimates <- paste(datafolder, "Harvard_elbow2_allyears_all.Rdata", sep = "")
save.image(file = file_estimates )








