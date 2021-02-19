# examples with multiple covariates for paper (version October 2020)

rm(list=ls())


require(mclust)
require(RSpectra)
require(grdpg)


mc_simulations <- function(index, n, Beta){
  #### Generate network
  ## Hyperparameter
  seed <- index
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
  #n <- 2000                  # Number of nodes
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
  #beta <- c(.5, .75)     # If two covariates: c(0,1, 0.3), etc
  beta <- Beta
  cov <- c(2,2)        # Possible value that the covariate could take, e.g. two binary = c(2,2)
  
  # Balanced case
  covariates <- matrix(NA,nrow = n, ncol = 2)
  covariates[,1] <- c(rep(1,n/4),rep(2,n/4), rep(1,n/4), rep(2,n/4))
  covariates[,2] <- c(rep(1,n/8),rep(2,n/4), rep(1,n/8),  rep(1,n/8),rep(2,n/4), rep(1,n/8))

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
  dhat <- 8
  
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
  
  result1 <- estimatebeta(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, sd=FALSE, link = 'logit', check = check)
  betahat1 <- sapply(result1$betahats, mean)

  result2 <- estimatebeta2(Xhat, muhats, Ipq, cov, covariates, clusters_cov, sd=FALSE, link = 'logit', check = check)
  betahat2 <- sapply(Map('*',result2$betahats,result2$pis), sum)

  runtime_GRDPG <- proc.time() - ptm
  filename <- paste("estimate",index,"_n", n, ".RData",sep = "")
  save.image(file = filename)
  
  output <- list(betahat1,betahat2,beta,runtime_GRDPG)
  names(output) <- c("betahat1","betahat2","beta","runtime_GRDPG")
  return(output)

  }



#### monte carlo
nsim <- 1000
n <- 10000
Beta <- c(0.5,0.75)
setwd("C:/Users/amele1/Dropbox/grdpg/paper/jasa_template/supplement/mc")
library(parallel)
number_cores <- 10
cl <- makeCluster(number_cores)
clusterEvalQ(cl, "library(grdpg)")
clusterEvalQ(cl, "require(RSpectra)")
clusterEvalQ(cl, "require(irlba)")
clusterEvalQ(cl, "require(mclust)")
clusterExport(cl, c("n", "Beta", "mc_simulations", "generateB", "generateP", "generateA", "SpectralEmbedding", "Mclust", "mclustBIC", "getClusters", "getIpq"))
clusterExport(cl, c("BFcheck", "Removecheck", "getBlockCovariates", "estimatebeta", "estimatebeta2", "logit", "sigmoid"))
results <- parLapply(cl, X=1:nsim, function(X) mc_simulations(X, n = n, Beta = Beta))
stopCluster(cl)


# make tables
estimates_grdpg <- data.frame(matrix(NA, nrow = nsim, ncol = 5))
names(estimates_grdpg) <- c("betahat1", "betahat2", "betahat1", "betahat2", "time")
for (table in 1:nsim){
  estimates_grdpg[table,1:2] <- results[[table]]$betahat1
  estimates_grdpg[table,3:4] <- results[[table]]$betahat2
  estimates_grdpg[table,5] <- results[[table]]$runtime_GRDPG[1]
}

# save session info
info_mc <- sessionInfo()

setwd("C:/Users/amele1/Dropbox/grdpg/paper/jasa_template/supplement/mc")
filename <- paste("mc_estimates_nsim", nsim, "_n", n, "_all.RData", sep="")
save.image(file = filename)



