#' @import Rcpp
#' @import abind
#' @import ggplot2
#' @import boot
#' @import bootstrap
#' @import DAAG
#' @import microbenchmark
#' @importFrom stats pchisq
#' @importFrom utils tail
#' @useDynLib BA23001021
NULL

#' @title A dataset used for illustration.
#' @name simdata
#' @description This dataset is simulated, with 2 time change-points at 3, 6
#' @examples
#' \dontrun{
#' data(simdata)
#' t <- simdata$t; delta <- simdata$delta; Z <- simdata$z
#' plik.est(t, delta, Z)
#' }
NULL

#' @title Log partial-likelihood.
#' @description Compute log partial-likelihood under Cox model, assume right censoring and no tie
#' @param t Survival time
#' @param delta Censoring indicator (1 failure, 0 right censoring)
#' @param Z Covariates
#' @param beta Coefficients
#' @return Log partial-likelihood (numeric)
#' @examples
#' \dontrun{
#' t <- rexp(1e3); delta <- rep(1, 1e3)
#' Z <- rnorm(1e3); beta <- 1
#' plik(t, delta, Z, beta)
#' }
#' @export
plik <- function(t, delta, Z, beta){
  Z <- cbind(Z); beta <- as.numeric(beta)
  ord <- order(t)
  t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]
  
  ss <- as.numeric(Z %*% beta)
  ee <- exp(ss)
  ff <- log(rev(cumsum(rev(ee))))
  
  return(sum(delta * (ss - ff)))
}

#' @title Derivative of log partial-likelihood.
#' @description Compute derivative of log partial-likelihood under Cox model, assume right censoring and no tie
#' @param t Survival time
#' @param delta Censoring indicator (1 failure, 0 right censoring)
#' @param Z Covariates
#' @param beta Coefficients
#' @return Derivative of log partial-likelihood (numeric vector)
#' @examples
#' \dontrun{
#' t <- rexp(1e3); delta <- rep(1, 1e3)
#' Z <- rnorm(1e3); beta <- 1
#' plik.1d(t, delta, Z, beta)
#' }
#' @export
plik.1d <- function(t, delta, Z, beta){
Z <- cbind(Z); beta <- as.numeric(beta)
ord <- order(t)
t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]

ss <- as.numeric(Z %*% beta)
ee <- exp(ss); ff <- rev(cumsum(rev(ee)))
return(p1d(ee, ff, delta, Z))
}

#' @title Second derivative of log partial-likelihood.
#' @description Compute second derivative of log partial-likelihood under Cox model, assume right censoring and no tie
#' @param t Survival time
#' @param delta Censoring indicator (1 failure, 0 right censoring)
#' @param Z Covariates
#' @param beta Coefficients
#' @return Second derivative of log partial-likelihood (numeric matrix)
#' @examples
#' \dontrun{
#' t <- rexp(1e3); delta <- rep(1, 1e3)
#' Z <- rnorm(1e3); beta <- 1
#' plik.2d(t, delta, Z, beta)
#' }
#' @export
plik.2d <- function(t, delta, Z, beta){
  Z <- cbind(Z); beta <- as.numeric(beta)
  ord <- order(t)
  t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]
  
  ss <- as.numeric(Z %*% beta)
  ee <- exp(ss); ff <- rev(cumsum(rev(ee)))
  
  return(p2d(ee, ff, delta, Z))
}

#' @title Partial-likelihood estimation.
#' @description Use Newton-Raphson method to compute maximum partial-likelihood estimation, assume right censoring and no tie
#' @param t Survival time
#' @param delta Censoring indicator (1 failure, 0 right censoring)
#' @param Z Covariates
#' @param init Initial value
#' @param control Convergence threshold
#' @return A list containing the estimation and the variance-covariance matrix
#' @examples
#' \dontrun{
#' t <- rexp(1e3); delta <- rep(1, 1e3)
#' Z <- rnorm(1e3)
#' plik.est(t, delta, Z)
#' }
#' @export
plik.est <- function(t, delta, Z, init = rep(0, NCOL(Z)), control = 1e-5){
  x <- init
  fx <- plik.1d(t, delta, Z, x)
  var <- chol2inv(chol(-plik.2d(t, delta, Z, x)))
  while(sum(abs(fx)) > control){
    x <- x + var %*% fx
    fx <- plik.1d(t, delta, Z, x)
    var <- chol2inv(chol(-plik.2d(t, delta, Z, x)))
  }
  
  return(list(est = as.numeric(x), var = var))
}

#' @title Detect multiple change-points by dynamic programming.
#' @description Use dynamic programming to detect multiple change-points
#' @param t Survival time
#' @param delta Censoring indicator (1 failure, 0 right censoring)
#' @param Z Covariates
#' @param candid Candid time points
#' @param ac 1 - significant level
#' @return A list containing all information
#' @export
dp.mcp <- function(t, delta, Z, candid, ac){
  len <- length(candid) - 1
  opt <- vector("list", len)
  
  for(i in len:1){
    l <- i; r <- len + 1
    
    delta.temp <- delta; delta.temp[t <= candid[l] | t > candid[r]] <- 0
    fit <- plik.est(t, delta.temp, Z); lik <- plik(t, delta.temp, Z, fit$est)
    beta <- rbind(fit$est)
    var <- fit$var; var <- array(var, dim = c(dim(var), 1))
    
    opt[[i]] <- list(
      index = 0, r = r,
      beta = beta, var = var, lik = lik
    )
    
    if(i < len) for(j in len:(i + 1)){
      l <- i; r <- j
      
      delta.temp <- delta; delta.temp[t <= candid[l] | t > candid[r]] <- 0
      fit <- plik.est(t, delta.temp, Z); lik <- plik(t, delta.temp, Z, fit$est)
      
      for (k in 1:length(opt[[j]]$lik)){
        rr <- opt[[j]]$r[k]
        
        beta2 <- opt[[j]]$beta[k, ]; var2 <- opt[[j]]$var[, , k]
        m <- fit$est - beta2; var <- fit$var + var2
        
        chistat <- as.numeric(m %*% chol2inv(chol(var)) %*% m)
        pvalue <- 1 - pchisq(q = chistat, df = NCOL(Z))
        
        if(pvalue < ac){
          lik <- lik + opt[[j]]$lik[k]
          opt[[i]]$index <- c(opt[[i]]$index, k)
          opt[[i]]$r <- c(opt[[i]]$r, r)
          opt[[i]]$beta <- rbind(opt[[i]]$beta, fit$est)
          opt[[i]]$var <- abind(opt[[i]]$var, fit$var, along = 3)
          opt[[i]]$lik <- c(opt[[i]]$lik, lik)
          break
        }
      }
    }
    o <- order(opt[[i]]$lik, decreasing = TRUE)
    opt[[i]]$index <- opt[[i]]$index[o]
    opt[[i]]$r <- opt[[i]]$r[o]
    opt[[i]]$beta <- opt[[i]]$beta[o, , drop = FALSE]
    opt[[i]]$var <- opt[[i]]$var[, , o, drop = FALSE]
    opt[[i]]$lik <- opt[[i]]$lik[o]
  }
  
  partition <- 1
  r <- opt[[1]]$r[1]
  index <- opt[[1]]$index[1]
  beta.res <- opt[[1]]$beta[1, , drop = FALSE]
  var.res <- opt[[1]]$var[, , 1, drop = FALSE]
  while ((r != len + 1) && (index != 0)) {
    partition <- c(partition, r)
    index0 <- opt[[r]]$index[1]
    r0 <- opt[[r]]$r[index]
    beta.res <- rbind(beta.res, opt[[r]]$beta[index, , drop = FALSE])
    var.res <- abind(var.res, opt[[r]]$var[, , index, drop = FALSE])
    r <- r0
    index <- index0
  }
  partition <- candid[partition]
  partition <- c(partition, tail(candid, 1))
  
  return(list(
    partition = partition, opt = opt, beta = beta.res, var = var.res
  ))
}
