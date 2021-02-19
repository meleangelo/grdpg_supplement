# examples for paper (version July 2019)

require(dplyr)
require(mclust)
require(RSpectra)
require(blockmodels)
require(ggplot2)
require(grdpg)


####################################################
########
######## Simulation Example 1 - without covariates
########
####################################################

# we generate networks with n=2000,5000,10000, 20000, 50000 nodes
# and K = 2,5,10 blocks, and d=1,2  

#######################################################
#### n=2000, K=2, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 2000                 # Number of nodes
K <- 2                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.7), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)








#######################################################
#### n=5000, K=2, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 5000                 # Number of nodes
K <- 2                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.7), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)




#######################################################
#### n=10000, K=2, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 10000                 # Number of nodes
K <- 2                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.7), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)




#######################################################
#### n=20000, K=2, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 20000                 # Number of nodes
K <- 2                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.7), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)





#######################################################
#### n=50000, K=2, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 50000                 # Number of nodes
K <- 2                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.7), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)





###########################
# K = 5
######




#######################################################
#### n=2000, K=5, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 2000                 # Number of nodes
K <- 5                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.3, 0.5, 0.7, 0.9), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)








#######################################################
#### n=5000, K=5, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 5000                 # Number of nodes
K <- 5                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.3, 0.5, 0.7, 0.9), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)




#######################################################
#### n=10000, K=5, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 10000                 # Number of nodes
K <- 5                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.3, 0.5, 0.7, 0.9), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)






#######################################################
#### n=20000, K=5, d=1
#######################################################

#### Generate network
## Hyperparameter
seed <- 2018
addCovariates <- FALSE
dmax <- 5                # Maximal embedding dimension
n <- 20000                 # Number of nodes
K <- 5                    # Number of blocks
d <- 1                    # Dimension of latent position

## Latent position
latent <- matrix(cbind(0.1, 0.3, 0.5, 0.7, 0.9), nrow = d)                         # d = 1


## Balanced case
pi <- rep(1/K, K)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}


## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, P, directed = FALSE, type = "bernoulli", seed)


## Estimation with GRDPG approach
ptm <- proc.time()

## ASE
set.seed(1977)
embed <- SpectralEmbedding(A, dmax, work = 100)
s <- embed$D

## Choose embed dimension
dimselect(s)
dhat <- dimselect(s)$elbow[1]

## Estimate latent position
Xhat <- embed$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))
#temp <- eigen(A)
Ipq <- getIpq(A, dhat)

## Cluster based on estimated latent position
set.seed(1977)
model <- Mclust(Xhat, G = 1:10)
clusters <- getClusters(data.frame(model$z))
muhats <- model$parameters$mean
muhats <- matrix(muhats, nrow = dhat)
# BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat <- t(muhats) %*% Ipq %*% muhats
BXhat
xxx <- SpectralEmbedding(BXhat, dhat, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)
runtime_GRDPG <- proc.time() - ptm


## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(est_latent)

## Visualization
# Screeplot
dat <- data.frame(s)
pp1 <- ggplot(dat, aes(x=1:dmax, y=s)) + geom_line() + geom_point() + scale_x_continuous(breaks = 1:dmax)
pp1 <- pp1 + labs(title = 'Screeplot (without Covariates)', x = 'Rank', y = 'Singular Value')
print(pp1)

# Latent position
dat <- data.frame(rotate(Xhat, latent, K))
if(d == 1) {
  names(dat) <- 'Xhat'
  Blocks <- as.factor(blocks)
  pp4 <- ggplot(dat, aes(Xhat)) + geom_histogram(aes(fill=Blocks), bins = 100)
  for (k in 1:length(latent)) {
    pp4 <- pp4 + geom_vline(xintercept = latent[k], color = 'black')
  }
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = '')
  print(pp4)
} else {
  Blocks <- as.factor(blocks)
  latent_vecs <- data.frame(t(latent))
  pp4 <- ggplot(dat) + geom_point(aes(X1, X2, color = Blocks), size = 0.3, alpha = 0.5)
  pp4 <- pp4 + geom_point(data = latent_vecs, aes(x=X1, y=X2), shape = 4, size = 5)
  pp4 <- pp4 + labs(title = 'Latent Position (without Covariates)', x = 'PC1', y = 'PC2')
  print(pp4)
}

multiplot(pp1, pp4, cols = 2)


# Estimation with VEM
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', A)
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxx <- SpectralEmbedding(Bhat, dhat, work = 100)
est_latentEM <- xxx$X * sqrt(xxx$D)
runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(est_latentEM)












########################################################
###
### Example 3 - logit link
###
#######################################################

# n=2000,K=2, d=1, beta = 1.5

#######################################################
############## without covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- FALSE 
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 2000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, sigmoid(P), seed = seed)


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1], dhat)

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}
xxx <- SpectralEmbedding(BXhat, 1, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)

runtime_GRDPG <- proc.time() - ptm

## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(latent)        # estimated latent positions
print(est_latent)        # estimated latent positions

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
dat <- data.frame(Xhat, Blocks = as.factor(blocks))
pp4 <- ggplot(dat) + 
  geom_point(aes(X1, X2, color = Blocks), alpha = 0.5) + 
  labs(title = "Latent Position (without Covariates)", x = "PC1", y = "PC2")
multiplot(pp1, pp4, cols = 2)


#### Variational EM approach
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', as.matrix(A))
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxxEM <- SpectralEmbedding(logit(Bhat), 1, work = 100)
est_latentEM <- xxxEM$X * sqrt(xxxEM$D)

runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(logit(Bhat))    # Estimated B matrix
print(latent)
print(est_latentEM)



#######################################################
############## with covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 2000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
Icov <- matrix(0, nrow = n, ncol = n)
for (i in 1:nrow(Icov)) {
  for (j in 1:ncol(Icov)) {
    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
  }
}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[1]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)
multiplot(pp1, pp3, pp4, pp2[[1]], pp2[[2]], pp2[[3]], cols = 2)

#### Variational EM approach
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli_covariates('SBM_sym', as.matrix(A), list(Icov), explore_max = K*prod(cov)+1)
modelEM$estimate()

## Get clusters, estimated beta and B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
betahat <- modelEM$model_parameters[[Khat]]$beta
Bhat <- modelEM$model_parameters[[Khat]]$m
xxxEM <- SpectralEmbedding(logit(Bhat), 1, work = 100)
est_latentEM <- xxxEM$X * sqrt(xxxEM$D)

runtime_EM <- proc.time() - ptm

# For K = 2 with one binary covariate
# clustersEM2 <- getClusters(data.frame(modelEM$memberships[[2]]$Z))
# adjustedRandIndex(blocks, clustersEM2)
# clustersEM4 <- getClusters(data.frame(modelEM$memberships[[4]]$Z))
# adjustedRandIndex(blocks_cov, clustersEM4)
# Bhat4 <- modelEM$model_parameters[[4]]$m

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(beta)           # True beta
print(betahat)        # Estimated beta
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(latent)
print(est_latentEM)




#########################################################






# n=5000,K=2, d=1, beta = 1.5

#######################################################
############## without covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- FALSE 
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 5000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, sigmoid(P), seed = seed)


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1], dhat)

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}
xxx <- SpectralEmbedding(BXhat, 1, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)

runtime_GRDPG <- proc.time() - ptm

## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(latent)        # estimated latent positions
print(est_latent)        # estimated latent positions

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
dat <- data.frame(Xhat, Blocks = as.factor(blocks))
pp4 <- ggplot(dat) + 
  geom_point(aes(X1, X2, color = Blocks), alpha = 0.5) + 
  labs(title = "Latent Position (without Covariates)", x = "PC1", y = "PC2")
multiplot(pp1, pp4, cols = 2)


#### Variational EM approach
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli('SBM_sym', as.matrix(A))
modelEM$estimate()

## Get clusters and estimated B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
Bhat <- modelEM$model_parameters[[Khat]]$pi
xxxEM <- SpectralEmbedding(logit(Bhat), 1, work = 100)
est_latentEM <- xxxEM$X * sqrt(xxxEM$D)

runtime_EM <- proc.time() - ptm

## Evaluation
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(B)              # True B matrix
print(logit(Bhat))    # Estimated B matrix
print(latent)
print(est_latentEM)



#######################################################
############## with covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 5000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
Icov <- matrix(0, nrow = n, ncol = n)
for (i in 1:nrow(Icov)) {
  for (j in 1:ncol(Icov)) {
    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
  }
}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[1]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)


#### Variational EM approach
ptm <- proc.time()

## Estimation
modelEM <- BM_bernoulli_covariates('SBM_sym', as.matrix(A), list(Icov), explore_max = K*prod(cov)+1)
modelEM$estimate()

## Get clusters, estimated beta and B matrix
Khat <- which.max(modelEM$ICL)
clustersEM <- getClusters(data.frame(modelEM$memberships[[Khat]]$Z))
betahat <- modelEM$model_parameters[[Khat]]$beta
Bhat <- modelEM$model_parameters[[Khat]]$m
xxxEM <- SpectralEmbedding(logit(Bhat), 1, work = 100)
est_latentEM <- xxxEM$X * sqrt(xxxEM$D)

runtime_EM <- proc.time() - ptm

# For K = 2 with one binary covariate
# clustersEM2 <- getClusters(data.frame(modelEM$memberships[[2]]$Z))
# adjustedRandIndex(blocks, clustersEM2)
# clustersEM4 <- getClusters(data.frame(modelEM$memberships[[4]]$Z))
# adjustedRandIndex(blocks_cov, clustersEM4)
# Bhat4 <- modelEM$model_parameters[[4]]$m

## Evaluation
sink("Users/Angelo/Dropbox/DNW/codes/EM_n5000_K2_withcov.txt")
ARI_EM <- adjustedRandIndex(blocks, clustersEM)
print(ARI_EM)         # ARI
print(runtime_EM)     # Runtime
print(beta)           # True beta
print(betahat)        # Estimated beta
print(B)              # True B matrix
print(Bhat)           # Estimated B matrix
print(latent)
print(est_latentEM)
sink()












#########################################################






# n=10000,K=2, d=1, beta = 1.5 ONLY GRDPG

#######################################################
############## without covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- FALSE 
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 10000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates)
P <- generateP(latent, d, block_size, addCovariates)
A <- generateA(n, sigmoid(P), seed = seed)


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- ifelse(is.null(dhat), dimselect(s)$elbow[1], dhat)

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}
xxx <- SpectralEmbedding(BXhat, 1, work = 100)
est_latent <- xxx$X * sqrt(xxx$D)

runtime_GRDPG <- proc.time() - ptm

## Evaluation
ARI_GRDPG <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG)         # ARI
print(runtime_GRDPG)     # Runtime
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix
print(latent)        # estimated latent positions
print(est_latent)        # estimated latent positions

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
dat <- data.frame(Xhat, Blocks = as.factor(blocks))
pp4 <- ggplot(dat) + 
  geom_point(aes(X1, X2, color = Blocks), alpha = 0.5) + 
  labs(title = "Latent Position (without Covariates)", x = "PC1", y = "PC2")
multiplot(pp1, pp4, cols = 2)




#######################################################
############## with covariates
######################################################

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 10000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 1                     # Dimension of latent position
latent <- cbind(-1.5, 1)   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[1]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)






###################################################
##
## EXAMPLE WITH d=2
##
##################################################



#######################################################
############## with covariates
######################################################

### ##### n = 2000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 2000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[2]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)




#######################################################
############## with covariates
######################################################

### ##### n = 5000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 5000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[2]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)







#######################################################
############## with covariates
######################################################

### ##### n = 10000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 10000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 1.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm

## Post analysis
Qhat <- Xhat %*% Ipq %*% t(Xhat)
if (check == 'BF') {
  What <- logit(BFcheck(Qhat))
} else {
  What <- logit(Removecheck(Qhat))
}
Aprime <- getAwithoutCovariates(What, betahat1, covariates)

embedprime <- SpectralEmbedding(Aprime, dmax, maxit = maxit, work = work, tol = tol)
sprime_simple <- embedprime$D
dhatprime <- dimselect(sprime_simple)$elbow[2]
Xhatprime <- embedprime$X[,1:dhatprime] %*% sqrt(diag(sprime_simple[1:dhatprime], nrow=dhatprime, ncol=dhatprime))

model2 <- Mclust(Xhatprime, G, verbose = FALSE)
clusters <- getClusters(data.frame(model2$z))




## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
ARI_GRDPG2 <- adjustedRandIndex(blocks, clusters)
print(ARI_GRDPG1)        # ARI with covariates
print(ARI_GRDPG2)        # ARI without covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(betahat2)          # Estimated beta
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix

## Visualization
cols <- ncol(A)
temp1 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(A), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp1 <- scree(tempdat$raw[1:dmax], 'Screeplot (with Covariates)')
pp2 <- plotLatentPosition(Xhat, blocks, TRUE, dhat, covariates)

cols <- ncol(Aprime)
temp1 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'LA')
s1 <- temp1$values
temp2 <- eigs_sym(matrix(as.numeric(Aprime), ncol = cols), dmax, 'SA')
s2 <- temp2$values
tempdat <- data.frame(raw = c(s1,s2)) %>%
  mutate(sign = ifelse(raw>0, 'positive', 'negative'), s = abs(raw)) %>%
  arrange(desc(s))
pp3 <- scree(tempdat$raw[1:dmax], 'Screeplot (without Covariates)')
pp4 <- plotLatentPosition(Xhatprime, blocks, FALSE, latent, K, d)
multiplot(pp1, pp3, pp2[[1]], pp4, cols = 2)






###################################################
##
## EXAMPLE WITH d=2 and small beta
##
##################################################



#######################################################
############## with covariates
######################################################

### ##### n = 2000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 2000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 0.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm



## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
print(ARI_GRDPG1)        # ARI with covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(sqrt(sd21)/n)      # std error
print(betahat2)          # Estimated beta
print(sqrt(sd22)/n)      # std error
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix



#######################################################
############## with covariates
######################################################

### ##### n = 5000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 5000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 0.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm



## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
print(ARI_GRDPG1)        # ARI with covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(sqrt(sd21)/n)      # std error
print(betahat2)          # Estimated beta
print(sqrt(sd22)/n)      # std error
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix







#######################################################
############## with covariates
######################################################

### ##### n = 10000

#### Generate network
## Hyperparameter
seed <- 201906
addCovariates <- TRUE
G <- 1:9                 # Used for GMM. If G = 1:9, then GMM will consider all the combination of 1 block to 9 block.
dmax <- 10               # Used for embedding. If dmax = 10, then will only calculate the first 10 singular value in SVD.
dhat <- NULL             # Used for estimating latent position. If dhat = NULL, then choose dhat using profile likelihood.
maxit <- 1000            # Used for embedding, specifically, `irlba` function for SVD. 
work <- 50               # Used for embedding, specifically, `irlba` function for SVD. 
tol <- 1e-05             # Used for embedding, specifically, `irlba` function for SVD. 
check <- 'BF'            # Used for checking probability

## Parameter
dmax <- 10                 # Maximal embedding dimension
n <- 10000                  # Number of nodes
K <- 2                     # Number of blocks
d <- 2                     # Dimension of latent position
latent <- cbind(c(-1.5,-1), c(1, .5))   # Latent position                   

## Balanced case
pi_1 <- 0.5
pi_2 <- 1 - pi_1
pi <- c(pi_1, pi_2)
block_size <- round(pi * n)
blocks <- c()
for (k in 1:length(block_size)) {
  blocks <- c(blocks, rep(k, block_size[k]))
}

## One binary covariate
# For (G)RDPG
beta <- 0.5     # If two covariates: c(0,1, 0.3), etc
cov <- 2        # Possible value that the covariate could take, e.g. two binary = c(2,2)

# Balanced case
pi_z_11 <- 0.5
pi_z_12 <- 1 - pi_z_11
pi_z_21 <- 0.5
pi_z_22 <- 1 - pi_z_21
pi_z <- c(pi_z_11, pi_z_12, pi_z_21, pi_z_22)
pi_cov <- c(pi[1]*pi_z[1], pi[1]*pi_z[2], pi[2]*pi_z[3], pi[2]*pi_z[4])  
block_size_cov <- round(pi_cov * n)
covariates <- matrix(c(rep(1,block_size_cov[1]), rep(2,block_size_cov[2]), rep(1,block_size_cov[3]), rep(2,block_size_cov[4])))
blocks_cov <- c()
for (k in 1:length(block_size_cov)) {
  blocks_cov <- c(blocks_cov, rep(k, block_size_cov[k]))
}

# For Variational EM
#Icov <- matrix(0, nrow = n, ncol = n)
#for (i in 1:nrow(Icov)) {
#  for (j in 1:ncol(Icov)) {
#    Icov[i,j] <- ifelse(covariates[i]==covariates[j], 1, 0)
#  }
#}

## Generate network (adjacency matrix) from (G)RDPG
B <- generateB(latent, K, d, addCovariates, cov, beta)                      
P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
A <- generateA(n, sigmoid(P), seed = seed)      


#### (G)RDPG approach
ptm <- proc.time()

## Embedding (giving adjacency matrix A)
embedding <- SpectralEmbedding(A, dmax, maxit = maxit, work = work, tol = tol)

## Choose embed dimension, i.e. dhat, using profile likelihood
s <- embedding$D
dhat <- 4

## Construct I_pq matrix for GRDPG (do not need this for RDPG)
Ipq <- getIpq(A, dhat)

## Estimate latent position
Xhat <- embedding$X[,1:dhat] %*% sqrt(diag(s[1:dhat], nrow=dhat, ncol=dhat))

## Cluster based on estimated latent position
model <- Mclust(Xhat, G, verbose = FALSE)
clusters_cov <- getClusters(data.frame(model$z))

## Estimate block probability matrix
muhats <- matrix(model$parameters$mean, nrow = dhat)
BXhat <- t(muhats) %*% Ipq %*% muhats
if (check == 'BF') {
  BXhat <- logit(BFcheck(BXhat))
} else {
  BXhat <- logit(Removecheck(BXhat))
}

## Estimate beta
covariates_block <- getBlockCovariates(covariates, clusters_cov)

result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = check)
betahat1 <- sapply(result1$betahats, mean)
betahat1_unbiased <- betahat1 - sapply(result1$bias, mean)
sd21 <- sapply(result1$sd2s, mean)

result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'logit', check = check)
betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)
betahat2_unbiased <- betahat2 - sapply(Map('*',result2$bias,result2$pis), sum)
sd22 <- sapply(Map('*',result2$sd2s,result2$pis), sum)

runtime_GRDPG <- proc.time() - ptm



## Evaluation
ARI_GRDPG1 <- adjustedRandIndex(blocks_cov, clusters_cov)
print(ARI_GRDPG1)        # ARI with covariates
print(runtime_GRDPG)     # Runtime
print(beta)              # True beta
print(betahat1)          # Estimated beta
print(sqrt(sd21)/n)      # std error
print(betahat2)          # Estimated beta
print(sqrt(sd22)/n)      # std error
print(B)                 # True B matrix
print(BXhat)             # Estimated B matrix


