### TESTING VARIANCE FORMULAS WITH PLUGIN ESTIMATORS

# examples with multiple covariates for paper (version October 2020)
# this simulation assumes that the unobserved blocks are unbalanced;
# the two binary covariates are independent and balanced

rm(list=ls())


require(mclust)
require(RSpectra)
require(grdpg)
require(dplyr)


#### monte carlo
nsim <- 100
n <- 10000
Beta <- c(0.5,0.75)
setwd("/Users/Angelo/Dropbox/grdpg/paper/jasa_template/supplement/mc")
library(parallel)
number_cores <- 5
cl <- makeCluster(number_cores)
clusterEvalQ(cl, "require(grdpg)")
clusterEvalQ(cl, "require(RSpectra)")
clusterEvalQ(cl, "require(irlba)")
clusterEvalQ(cl, "require(mclust)")
clusterEvalQ(cl, "require(dplyr)")
clusterEvalQ(cl, "n()")
require(mclust)
require(RSpectra)
require(grdpg)
require(dplyr)


getWeight <- function (clusters) 
{
  n <- length(clusters)
  weight <- data.frame(clusters) %>% group_by(clusters) %>% 
    summarise(count = n()) %>% mutate(freq = count/n)
  return(weight)
}
compute_cov <- function(i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq) {
  .Call("_grdpg_compute_cov", PACKAGE = "grdpg", i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq)
}

estimatebeta_T <- function(Xhat, muhats, Ipq, cov, covariates_block, clusters_cov, link = 'logit', check = 'BF', sd = TRUE, rho = 1) {
  require(mclust)
  require(RSpectra)
  require(grdpg)
  require(dplyr)
  result <- list()
  
  covariates_block <- data.frame(covariates_block)
  if (length(cov) != ncol(covariates_block)) {
    stop("The length of `cov` should equal to the number of columns in `covariates_block` (both equal to the number of covariates).")
  }
  if (!(link %in% c('identity', 'logit'))) {
    print("Unrecognized `link`, would use 'identity' by default.")
  }
  if (!(check %in% c('BF', 'Remove'))) {
    print("Unrecognized `check`, would use 'BF' by default.")
  }
  
  K <- length(unique(clusters_cov))
  BXhat <- t(muhats) %*% Ipq %*% muhats
  if (link == 'logit') {
    if (check == 'BF') {
      BXhat <- logit(BFcheck(BXhat))
    } else {
      BXhat <- logit(Removecheck(BXhat))
    }
  }
  betahats <- vector('list', ncol(covariates_block))
  if (sd) {
    theta <- t(muhats) %*% Ipq %*% muhats
    theta <- BFcheck(theta)
    eta <- getWeight(clusters_cov)$freq
    delta <- 0
    for (kk in 1:length(eta)) {
      delta <- delta + eta[kk] * muhats[,kk,drop=FALSE] %*% t(muhats[,kk,drop=FALSE])
    }
    deltainv <- solve(delta)
    zeta <- t(muhats) %*% deltainv %*% muhats
    E <- vector('list', length(unique(clusters_cov)))
    for (alpha in unique(clusters_cov)) {
      for (n in 1:nrow(Xhat)) {
        if (rho == 0) {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        } else {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(1-theta[alpha,clusters_cov[n]])*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        }
      }
    }
    bias <- vector('list', ncol(covariates_block))
    sd2s <- vector('list', ncol(covariates_block))
    psiil1s <- vector('list', ncol(covariates_block))
    psiil2s <- vector('list', ncol(covariates_block))
    sigma2il1s <- vector('list', ncol(covariates_block))
    sigma2il2s <- vector('list', ncol(covariates_block))
    sigma2l1l1s <- vector('list', ncol(covariates_block))
    sigma2l2l2s <- vector('list', ncol(covariates_block))
    covil1il2s <- vector('list', ncol(covariates_block))
  }
  #model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)
  #c <- getClusters(data.frame(model2$z))
  #c <- kmeans(diag(BXhat), ncol(BXhat)/prod(cov))$cluster
  c <- c(1,1,1,1,2,2,2,2)
  for (i in 1:nrow(BXhat)) {
    for (k in 1:ncol(covariates_block)) {
      ind1 <- which(covariates_block[,k]==covariates_block[i,k])
      ind2 <- which(covariates_block[,k]!=covariates_block[i,k])
      for (l1 in ind1) {
        for (l2 in ind2) {
          if (c[l1] == c[l2]) {
            temp <- setdiff(1:ncol(covariates_block), k)
            ind <- c()
            if (length(temp) > 0) {
              for (l in temp) {
                ind <- c(ind, covariates_block[l1,l] == covariates_block[l2,l])
              }
            } else {
              ind <- TRUE
            }
            if (all(ind)) {
              betahats[[k]] <- c(betahats[[k]], BXhat[i,l1] - BXhat[i,l2])
              if (sd) {
                meanE <- sapply(E, mean)
                covs_temps <- compute_cov(i, l1, l2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq)
                covs <- covs_temps[[1]]
                temps <- covs_temps[[2]]
                psiil1 <- temps[1]
                psiil2 <- temps[2]
                sigma2il1 <- temps[3]
                sigma2il2 <- temps[4]
                sigma2l1l1 <- temps[5]
                sigma2l2l2 <- temps[6]
                covil1il2 <- covs[1] + covs[2] + covs[3] + covs[4] - covs[5] - covs[6] - covs[7] - covs[8] + covs[9]
                if (link == 'logit') {
                  psiil1 <- psiil1 * logitderivative(BFcheck(theta[i,l1]))
                  psiil2 <- psiil2 * logitderivative(BFcheck(theta[i,l2]))
                  sigma2il1 <- sigma2il1 * logitderivative(BFcheck(theta[i,l1]))^2
                  sigma2il2 <- sigma2il2 * logitderivative(BFcheck(theta[i,l2]))^2
                  sigma2l1l1 <- sigma2l1l1 * logitderivative(BFcheck(theta[l1,l1]))^2
                  sigma2l2l2 <- sigma2l2l2 * logitderivative(BFcheck(theta[l2,l2]))^2
                  covil1il2 <- covil1il2 * logitderivative(BFcheck(theta[i,l1])) * logitderivative(BFcheck(theta[i,l2]))
                }
                psiil1s[[k]] <- c(psiil1s[[k]], psiil1)
                psiil2s[[k]] <- c(psiil1s[[k]], psiil2)
                sigma2il1s[[k]] <- c(sigma2il1s[[k]], sigma2il1)
                sigma2il2s[[k]] <- c(sigma2il2s[[k]], sigma2il2)
                sigma2l1l1s[[k]] <- c(sigma2l1l1s[[k]], sigma2l1l1)
                sigma2l2l2s[[k]] <- c(sigma2l2l2s[[k]], sigma2l2l2)
                covil1il2s[[k]] <- c(covil1il2s[[k]], covil1il2)
                if (i == l1) {
                  tempsd <- sigma2l1l1 + sigma2il2 - 2 * covil1il2
                } else if (i == l2) {
                  tempsd <- sigma2il1 + sigma2l2l2 - 2 * covil1il2
                } else {
                  tempsd <- sigma2il1 + sigma2il2 - 2 * covil1il2
                }
                tempbias <- (psiil1 - psiil2) / length(clusters_cov)
                bias[[k]] <- c(bias[[k]], tempbias)
                sd2s[[k]] <- c(sd2s[[k]], tempsd)
              }
            }
          }
        }
      }
    }
  }
  result$betahats <- betahats
  if (sd) {
    result$bias <- bias
    result$sd2s <- sd2s
    result$psiil1s <- psiil1s
    result$psiil2s <- psiil2s
    result$sigma2il1s <- sigma2il1s
    result$sigma2il2s <- sigma2il2s
    result$sigma2l1l1s <- sigma2l1l1s
    result$sigma2l2l2s <- sigma2l2l2s
    result$covil1il2s <- covil1il2s
  }
  
  return(result)
  
}

estimatebeta_TW <- function(Xhat, muhats, Ipq, cov, covariates, clusters_cov, link = 'identity', check = 'BF', sd = TRUE, rho = 1) {
  require(mclust)
  require(RSpectra)
  require(grdpg)
  require(dplyr)
  result <- list()
  
  K <- length(unique(clusters_cov))
  BXhat <- t(muhats) %*% Ipq %*% muhats
  if (link == 'logit') {
    if (check == 'BF') {
      BXhat <- logit(BFcheck(BXhat))
    } else {
      BXhat <- logit(Removecheck(BXhat))
    }
  }
  betahats <- vector('list', ncol(covariates))
  pis <- vector('list', length(cov))
  if (sd) {
    theta <- t(muhats) %*% Ipq %*% muhats
    theta <- BFcheck(theta)
    eta <- getWeight(clusters_cov)$freq
    delta <- 0
    for (kk in 1:length(eta)) {
      delta <- delta + eta[kk] * muhats[,kk,drop=FALSE] %*% t(muhats[,kk,drop=FALSE])
    }
    deltainv <- solve(delta)
    zeta <- t(muhats) %*% deltainv %*% muhats
    E <- vector('list', length(unique(clusters_cov)))
    for (alpha in unique(clusters_cov)) {
      for (n in 1:nrow(Xhat)) {
        if (rho == 0) {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        } else {
          E[[alpha]] <- c(E[[alpha]], theta[alpha,clusters_cov[n]]*(1-theta[alpha,clusters_cov[n]])*(Xhat[n,,drop=FALSE]%*%t(Xhat[n,,drop=FALSE])))
        }
      }
    }
    bias <- vector('list', ncol(covariates))
    sd2s <- vector('list', ncol(covariates))
    psiij1s <- vector('list', ncol(covariates))
    psiij2s <- vector('list', ncol(covariates))
    sigma2ij1s <- vector('list', ncol(covariates))
    sigma2ij2s <- vector('list', ncol(covariates))
    sigma2j1j1s <- vector('list', ncol(covariates))
    sigma2j2j2s <- vector('list', ncol(covariates))
    covij1ij2s <- vector('list', ncol(covariates))
  }
  #model2 <- Mclust(diag(BXhat), ncol(BXhat)/prod(cov), verbose = FALSE)
  #c <- getClusters(data.frame(model2$z))
  #c <- kmeans(diag(BXhat), ncol(BXhat)/prod(cov))$cluster
  c <- c(1,1,1,1,2,2,2,2)
  weight_cov <- getWeight(clusters_cov)
  weight <- data.frame()
  for (clusters in unique(c)) {
    ind <- which(c==clusters)
    count <- sum(weight_cov$count[weight_cov$clusters %in% ind])
    freq <- sum(weight_cov$freq[weight_cov$clusters %in% ind])
    weight <- rbind(weight, data.frame(clusters,count,freq))
  }
  weights_cov <- vector('list', ncol(covariates))
  for (k in 1:ncol(covariates)) {
    temp <- data.frame(clusters_cov, covariates[,k])
    names(temp)[2] <- 'covariates'
    weights_cov[[k]] <- temp %>%
      group_by(clusters_cov, covariates) %>%
      summarise(count = n()) %>%
      left_join(weight_cov, by = c('clusters_cov'='clusters')) %>%
      mutate(freq = count.x/count.y)
  }
  for (k in 1:ncol(covariates)) {
    for (i in 1:nrow(BXhat)) {
      for (j1 in 1:(ncol(BXhat)-1)) {
        for (j2 in (j1+1):ncol(BXhat)) {
          if (c[j1] == c[j2]) {
            t1 <- filter(weights_cov[[k]], clusters_cov == j1)
            t2 <- filter(weights_cov[[k]], clusters_cov == j2)
            for (l1 in 1:nrow(t1)) {
              ind1 <- setdiff(unique(t2$covariates), t1$covariates[l1])
              if (length(ind1) > 0) {
                for (l2 in ind1) {
                  # pi_1 <- weight$freq[weight$clusters==c[i]]
                  # pi_2 <- weight$freq[weight$clusters==c[j1]]
                  pi_cov_1 <- weight_cov$freq[weight_cov$clusters==j1]
                  pi_cov_2 <- weight_cov$freq[weight_cov$clusters==j2]
                  pi_z_1 <- t1$freq[l1]
                  pi_z_2 <- t2$freq[t2$covariates==l2]
                  temp <- setdiff(1:ncol(covariates), k)
                  pi_0 <- 1
                  if (length(temp) > 0) {
                    for (l in temp) {
                      temp_pi_0 <- 0
                      t3 <- filter(weights_cov[[l]], clusters_cov == j1)
                      t4 <- filter(weights_cov[[l]], clusters_cov == j2)
                      for (t_cov in unique(t3$covariates)) {
                        if (t_cov %in% t4$covariates) {
                          temp_pi_0 <- temp_pi_0 + sum(t3$freq[t3$covariates==t_cov]) * sum(t4$freq[t4$covariates==t_cov])
                        } else {
                          temp_pi_0 <- temp_pi_0 + 0
                        }
                      }
                      pi_0 <- c(pi_0, temp_pi_0)
                    }
                  }
                  betahats[[k]] <- c(betahats[[k]], BXhat[i,j1] - BXhat[i,j2])
                  pis[[k]] <- c(pis[[k]], prod(pi_0) * pi_cov_1 * pi_cov_2 * pi_z_1 * pi_z_2)
                  if (sd) {
                    meanE <- sapply(E, mean)
                    covs_temps <- compute_cov(i, j1, j2, K, meanE, muhats, deltainv, eta, theta, zeta, Ipq)
                    covs <- covs_temps[[1]]
                    temps <- covs_temps[[2]]
                    psiij1 <- temps[1]
                    psiij2 <- temps[2]
                    sigma2ij1 <- temps[3]
                    sigma2ij2 <- temps[4]
                    sigma2j1j1 <- temps[5]
                    sigma2j2j2 <- temps[6]
                    covij1ij2 <- covs[1] + covs[2] + covs[3] + covs[4] - covs[5] - covs[6] - covs[7] - covs[8] + covs[9]
                    if (link == 'logit') {
                      psiij1 <- psiij1 * logitderivative(BFcheck(theta[i,j1]))
                      psiij2 <- psiij2 * logitderivative(BFcheck(theta[i,j2]))
                      sigma2ij1 <- sigma2ij1 * logitderivative(BFcheck(theta[i,j1]))^2
                      sigma2ij2 <- sigma2ij2 * logitderivative(BFcheck(theta[i,j2]))^2
                      sigma2j1j1 <- sigma2j1j1 * logitderivative(BFcheck(theta[j1,j1]))^2
                      sigma2j2j2 <- sigma2j2j2 * logitderivative(BFcheck(theta[j2,j2]))^2
                      covij1ij2 <- covij1ij2 * logitderivative(BFcheck(theta[i,j1])) * logitderivative(BFcheck(theta[i,j2]))
                    }
                    psiij1s[[k]] <- c(psiij1s[[k]], psiij1)
                    psiij2s[[k]] <- c(psiij1s[[k]], psiij2)
                    sigma2ij1s[[k]] <- c(sigma2ij1s[[k]], sigma2ij1)
                    sigma2ij2s[[k]] <- c(sigma2ij2s[[k]], sigma2ij2)
                    sigma2j1j1s[[k]] <- c(sigma2j1j1s[[k]], sigma2j1j1)
                    sigma2j2j2s[[k]] <- c(sigma2j2j2s[[k]], sigma2j2j2)
                    covij1ij2s[[k]] <- c(covij1ij2s[[k]], covij1ij2)
                    if (i == j1) {
                      tempsd <- sigma2j1j1 + sigma2ij2 - 2 * covij1ij2
                    } else if (i == j2) {
                      tempsd <- sigma2ij1 + sigma2j2j2 - 2 * covij1ij2
                    } else {
                      tempsd <- sigma2ij1 + sigma2ij2 - 2 * covij1ij2
                    }
                    tempbias <- (psiij1 - psiij2) / length(clusters_cov)
                    bias[[k]] <- c(bias[[k]], tempbias)
                    sd2s[[k]] <- c(sd2s[[k]], tempsd)
                  }
                }
              }
            }
          }
        }
      }
    }
    symbol <- ifelse(BXhat[1,1]-BXhat[1,which(c==c[1])[2]]>0, 1, -1)
    pis[[k]] <- pis[[k]] / sum(pis[[k]])
    betahats[[k]] <- abs(betahats[[k]]) * symbol
  }
  result$betahats <- betahats
  result$pis <- pis
  if (sd) {
    result$bias <- bias
    result$sd2s <- sd2s
    result$psiij1s <- psiij1s
    result$psiij2s <- psiij2s
    result$sigma2ij1s <- sigma2ij1s
    result$sigma2ij2s <- sigma2ij2s
    result$sigma2j1j1s <- sigma2j1j1s
    result$sigma2j2j2s <- sigma2j2j2s
    result$covij1ij2s <- covij1ij2s
  }
  
  return(result)
  
}



mc_simulations <- function(index, n, Beta){
  require(mclust)
  require(RSpectra)
  require(grdpg)
  require(dplyr)
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
  pi_1 <- 0.3
  pi_2 <- 1 - pi_1
  pi <- c(pi_1, pi_2)
  block_size <- round(pi * n)
  blocks <- c()
  for (k in 1:length(block_size)) {
    blocks <- c(blocks, rep(k, block_size[k]))
  }
  
  
  ## two binary covariate
  # For (G)RDPG
  #beta <- c(.5, .75)     # If two covariates: c(0,1, 0.3), etc
  beta <- Beta
  cov <- c(2,2)        # Possible value that the covariate could take, e.g. two binary = c(2,2)
  
  # Balanced case
  covariates <- matrix(NA,nrow = n, ncol = 2)
  covariates[,1] <- rbinom(n,1,.5)
  covariates[,2] <- rbinom(n,1,.5)
  
  ## Generate network (adjacency matrix) from (G)RDPG
  B <- generateB(latent, K, d, addCovariates, cov, beta)                      
  P <- generateP(latent, d, block_size, addCovariates, covariates, beta)   
  A <- generateA(n, sigmoid(P), seed = seed)      
  
  # true latent positions for the model with covariates
  # svdB <- svd(sigmoid(B))
  # muhats_true <- svdB$u %*% sqrt(diag(svdB$d, nrow=8, ncol=8))
  svdB <- svd(sigmoid(B))
  muhats_true <- svdB$u %*% sqrt(diag(svdB$d))
  
  cluster_cov_true <- rep(NA, n)
  cluster_cov_true[blocks == 1 & covariates[,1] == 1 & covariates[,2] == 1] <- 1
  cluster_cov_true[blocks == 1 & covariates[,1] == 1 & covariates[,2] == 0] <- 2
  cluster_cov_true[blocks == 1 & covariates[,1] == 0 & covariates[,2] == 1] <- 3
  cluster_cov_true[blocks == 1 & covariates[,1] == 0 & covariates[,2] == 0] <- 4
  cluster_cov_true[blocks == 2 & covariates[,1] == 1 & covariates[,2] == 1] <- 5
  cluster_cov_true[blocks == 2 & covariates[,1] == 1 & covariates[,2] == 0] <- 6
  cluster_cov_true[blocks == 2 & covariates[,1] == 0 & covariates[,2] == 1] <- 7
  cluster_cov_true[blocks == 2 & covariates[,1] == 0 & covariates[,2] == 0] <- 8
  
  Xtrue <- matrix(NA, nrow = n, ncol = 8)
  for (nu in 1:8) {
    Xtrue[cluster_cov_true==nu,] <- t(muhats_true[nu,])
  }
  
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
  covariates_block_true <- getBlockCovariates(covariates, cluster_cov_true)
  

  # true
  result_true <- estimatebeta_T(Xtrue, t(muhats_true), Ipq, cov, covariates_block_true, cluster_cov_true, check = "BF", link = "logit", sd = TRUE)
  betahat_true <- sapply(result_true$betahats, mean)
  betahat_true_unbiased <- betahat_true - sapply(result_true$bias, mean)
  var_true <- sapply(result_true$sd2s, mean)/dim(A)[1]
  sd_true <- sqrt(var_true)/dim(A)[1]
  print(betahat_true)
  betahat_true_unbiased
  sd_true
  
  result_trueWA <- estimatebeta_TW(Xtrue, t(muhats_true), Ipq, cov, covariates, cluster_cov_true,  link = "logit", sd = TRUE)
  betahat_trueWA <- sapply(Map('*',result_trueWA$betahats,result_trueWA$pis), sum)
  betahat_trueWA_unbiased <- betahat_trueWA - sapply(Map('*',result_trueWA$bias,result_trueWA$pis), sum)
  varweights_trueWA <- Map('*', result_trueWA$pis, result_trueWA$pis)
  var_trueWA <- sapply(Map('*',result_trueWA$sd2s,varweights_trueWA), sum)
  sd_trueWA <- sqrt(var_trueWA)/dim(A)[1]
  print(betahat_trueWA )
  betahat_trueWA_unbiased
  sd_trueWA
  
  
  result_WA <- estimatebeta_WA(Xhat, muhats, Ipq, cov, covariates, clusters_cov,  link = "logit", sd = TRUE)
  betahat_WA <- sapply(Map('*',result_WA$betahats,result_WA$pis), sum)
  betahat_WA_unbiased <- betahat_WA - sapply(Map('*',result_WA$bias,result_WA$pis), sum)
  varweights_WA <- Map('*', result_WA$pis, result_WA$pis)
  var_WA <- sapply(Map('*',result_WA$sd2s,varweights_WA), sum)
  sd_WA <- sqrt(var_WA)/dim(A)[1]
  print(betahat_WA)
  betahat_WA_unbiased
  sd_WA
  
  
  
  runtime_GRDPG <- proc.time() - ptm

  output <- list(result_true,result_trueWA, result_WA,beta,runtime_GRDPG)
  names(output) <- c("result_true", "result_trueWA","result_WA","beta","runtime_GRDPG")
  return(output)
  
}


clusterExport(cl, c("n", "Beta", "mc_simulations", "generateB", "generateP", "generateA", "SpectralEmbedding", "Mclust", "mclustBIC", "getClusters", "getIpq"))
clusterExport(cl, c("BFcheck", "Removecheck", "getBlockCovariates", "estimatebeta_WA", "estimatebeta_SA", "logit", "sigmoid"))
clusterExport(cl, c("estimatebeta_T", "%>%", "n", "group_by", "mutate", "summarise", "left_join", "estimatebeta_TW", "logitderivative", "getWeight", "compute_cov"))

results <- parLapply(cl, X=1:nsim, function(X) mc_simulations(X, n = n, Beta = Beta))
stopCluster(cl)


save.image(file = "VarianceMC_n10000.RData")

for (ii in 1:100){
  print(results[[ii]]$result_true$betahats[[1]][1])
}


# make table with all MC results (betas and std devs)
table <- matrix(NA, ncol = 12, nrow = nsim)

for (i in 1:nsim){
  result_true <- results[[i]]$result_true
  betahat_true <- sapply(result_true$betahats, mean)
  betahat_true_unbiased <- betahat_true - sapply(result_true$bias, mean)
  var_true <- sapply(result_true$sd2s, mean)/n
  sd_true <- sqrt(var_true)/n
  print(betahat_true)
  betahat_true_unbiased
  sd_true
  table[i,1:2] <- betahat_true
  table[i,7:8] <- sd_true
  
  result_trueWA <- results[[i]]$result_trueWA
  betahat_trueWA <- sapply(Map('*',result_trueWA$betahats,result_trueWA$pis), sum)
  betahat_trueWA_unbiased <- betahat_trueWA - sapply(Map('*',result_trueWA$bias,result_trueWA$pis), sum)
  varweights_trueWA <- Map('*', result_trueWA$pis, result_trueWA$pis)
  var_trueWA <- sapply(Map('*',result_trueWA$sd2s,varweights_trueWA), sum)
  sd_trueWA <- sqrt(var_trueWA)/n
  print(betahat_trueWA )
  betahat_trueWA_unbiased
  sd_trueWA
  table[i,3:4] <- betahat_trueWA
  table[i,9:10] <- sd_trueWA
  
  
  result_WA <- results[[i]]$result_WA
  betahat_WA <- sapply(Map('*',result_WA$betahats,result_WA$pis), sum)
  betahat_WA_unbiased <- betahat_WA - sapply(Map('*',result_WA$bias,result_WA$pis), sum)
  varweights_WA <- Map('*', result_WA$pis, result_WA$pis)
  var_WA <- sapply(Map('*',result_WA$sd2s,varweights_WA), sum)
  sd_WA <- sqrt(var_WA)/n
  print(betahat_WA)
  betahat_WA_unbiased
  sd_WA
  table[i,5:6] <- betahat_WA
  table[i,11:12] <- sd_WA
  
}

# name columns for visualization
table <- data.frame(table)
names(table) <- c("betahat1", "betahat2","betahat1w", "betahat2w","betahat1wa", "betahat2wa",
                  "sd1", "sd2","sd1w", "sd2w","sd1wa", "sd2wa")
# table with values for estimates exactly equal to beta (some have numerical error!)
table_true <- table[table[,1]>0.499999999,]

# make monte carlo table 
mc_table <- colMeans(table[,7:12])
mc_table_true <- colMeans(table_true[,7:12])
mc_table
mc_table_true
