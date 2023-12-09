## ----eval=FALSE---------------------------------------------------------------
#  plik <- function(t, delta, Z, beta){
#    Z <- cbind(Z); beta <- as.numeric(beta)
#    ord <- order(t)
#    t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]
#  
#    ss <- as.numeric(Z %*% beta)
#    ee <- exp(ss)
#    ff <- log(rev(cumsum(rev(ee))))
#  
#    return(sum(delta * (ss - ff)))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  plik.1d <- function(t, delta, Z, beta){
#  Z <- cbind(Z); beta <- as.numeric(beta)
#  ord <- order(t)
#  t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]
#  
#  ss <- as.numeric(Z %*% beta)
#  ee <- exp(ss); ff <- rev(cumsum(rev(ee)))
#  return(p1d(ee, ff, delta, Z))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  plik.2d <- function(t, delta, Z, beta){
#    Z <- cbind(Z); beta <- as.numeric(beta)
#    ord <- order(t)
#    t <- t[ord]; delta <- delta[ord]; Z <- Z[ord, , drop = FALSE]
#  
#    ss <- as.numeric(Z %*% beta)
#    ee <- exp(ss); ff <- rev(cumsum(rev(ee)))
#  
#    return(p2d(ee, ff, delta, Z))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  plik.est <- function(t, delta, Z, init = rep(0, NCOL(Z)), control = 1e-5){
#    x <- init
#    fx <- plik.1d(t, delta, Z, x)
#    var <- chol2inv(chol(-plik.2d(t, delta, Z, x)))
#    while(sum(abs(fx)) > control){
#      x <- x + var %*% fx
#      fx <- plik.1d(t, delta, Z, x)
#      var <- chol2inv(chol(-plik.2d(t, delta, Z, x)))
#    }
#  
#    return(list(est = as.numeric(x), var = var))
#  }

## -----------------------------------------------------------------------------
library(BA23001021)
data(simdata)
t <- simdata$t; delta <- simdata$delta; Z <- simdata$z

# 待选集为 t 的20等分点
candid <- c(-Inf, quantile(t, seq(0.05, 0.95, 0.05)), Inf)

# 使用 Bonferonni 来控制 FWER ， 1330 为 20 个点的动态规划所需的最大检验数
ac <- 0.05 / 1330

res <- dp.mcp(t, delta, Z, candid, ac)
res$partition

