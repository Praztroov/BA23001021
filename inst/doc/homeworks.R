## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
knitr::opts_chunk$set(fig.align = "center", fig.width = 7, fig.height=5)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 200
x <- rnorm(n, mean = 0, sd = 3)
y <- 2 + 3 * x + rnorm(n, mean = 0, sd = 2)
plot(x, y)

## -----------------------------------------------------------------------------
fit <- lm(y~x)
summary(fit)

## -----------------------------------------------------------------------------
plot(x, y)
abline(a = fit$coefficients[1], b = fit$coefficients[2])

## -----------------------------------------------------------------------------
knitr::kable(head(airquality), align = "c", caption = "测试表-airquality")

## -----------------------------------------------------------------------------
par(mfrow = c(1,2))
hist(airquality$Temp)
hist(airquality$Ozone)

## -----------------------------------------------------------------------------
library(ggplot2)
n <- 1000
x <- runif(n, -2 * pi, 2 * pi)
y <- sin(x) + rnorm(n, 0, 1)

p <- ggplot(
  data.frame(x = x, y = y),
  aes(x, y)) +
  geom_point() +
  stat_smooth(
    formula = y ~ x,
    method = loess, 
    aes(colour = "fitted"),
    linewidth = 1
    ) +
  stat_function(
    fun=sin,
    aes(colour = "origin"),
    linewidth = 1
    ) +
  guides(colour=guide_legend(title=NULL))
p

## -----------------------------------------------------------------------------
my.sample <- function(x, size, prob = NA){
  
  #默认概率均分
  if(length(prob) == 1) if(is.na(prob)){ 
    prob = rep(1, length(x)) / length(x)
  }
  
  if(
    sum(prob) != 1 || sum(prob < 0) > 0
  ) stop("invalid prob")
  
  df <- c(0, cumsum(prob))#分布函数
  u <- runif(size, 0, 1)
  location <- sapply(u, FUN = df.inv, df)
  
  return(x[location])
}

df.inv <- function(u, df){#分布函数的逆
  uu <- u < df
  a <- uu[-1]
  b <- (!uu)[-length(df)]
  
  return(which.max(a * b))
}


## -----------------------------------------------------------------------------
set.seed(123)

x <- c("a", "c", "f", "x", "z")
size <- 10000
prob <- c(0.1, 0.3, 0.2, 0.15, 0.25)

y <- my.sample(x, size, prob)

freq <- integer(length(x))

for(i in 1:5){
  freq[i] <- sum(y == x[i]) / size
}

print(freq - prob)

## ----echo=FALSE---------------------------------------------------------------
dd <- as.data.frame(matrix(NA, nrow = 2 * length(x), ncol = 3))

colnames(dd) <- c("group", "value", "x")

dd$group <- c(rep("prob", length(x)), rep("freq", length(x)))
dd$value <- c(prob, freq)
dd$x <- c(x, x)

p <- ggplot(
  data = dd,
  aes(x = x, y = value, fill = group)
  ) +
  geom_bar(
    stat = 'identity',
    colour = 'black',
    position='dodge'
  ) + 
  ggtitle("样本频率与原始频率对比图") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(fill = guide_legend(title=NULL)) +
  theme(legend.position = c(0.1, 0.7))
p

## -----------------------------------------------------------------------------
laplace <- function(size){
  u <- rbind(runif(size, min = 0, max = 1))
  
  x <- integer(length(u))
  uu <- u < 0.5
  x[uu] <- log(2 * u[uu])
  x[!uu] <- -log(2 * (1 - u[!uu]))
  
  return(x)
}

set.seed(123)
x <- laplace(1000)

## ----echo=FALSE---------------------------------------------------------------
dens <- function(x) return(0.5 * exp(-abs(x)))

p <- ggplot(data.frame(x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    colour = 'black',
    fill = "gray"
    ) +
  stat_function(
    fun = dens,
    linewidth = 1
  ) + 
  ggtitle("样本频率与目标密度对比图") +
  theme(plot.title = element_text(hjust = 0.5))
p

## ----echo=FALSE---------------------------------------------------------------
distr <- function(x) return(0.5 * exp(x) * (x < 0) + (1 - 0.5 * exp(-x)) * (x >= 0))

p <- ggplot(data.frame(x), aes(x = x)) +
  stat_function(
    fun = distr,
    aes(color = "target df"),
    xlim = c(min(x), max(x))
  )  +
  stat_ecdf(aes(color = "ecdf"))+ 
  ggtitle("样本累积经验分布与目标分布对比图") +
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(colour = guide_legend(title=NULL)) +
  theme(legend.position = c(0.8, 0.2))
p

## -----------------------------------------------------------------------------
beta.single <- function(a, b, c){
  while(TRUE){
    u <- runif(1, 0, 1)
    y <- runif(1, 0, 1)
    if(
      u <= y^(a-1) * (1-y)^(b-1) / c
    ) return(y)
  }
}

beta.sample <- function(size, a, b){
  if(a < 1 || b < 1) stop("invalid a or b")
  
  c <- (a-1)^(a-1) * (b-1)^(b-1) / (a+b-2)^(a+b-2)
  x <- sapply(rep(a, size), FUN = beta.single, b, c)
  
  return(x)
}

set.seed(123)
x <- beta.sample(1000, 3, 2)

## ----echo=FALSE---------------------------------------------------------------
dens <- function(x) return(x^2 * (1-x) / beta(3,2))

p <- ggplot(data.frame(x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    colour = 1,
    fill = "gray"
    ) +
  stat_function(
    fun = dens,
    linewidth = 1
  ) + 
  ggtitle("样本频率与理论密度对比图") +
  theme(plot.title = element_text(hjust = 0.5))
p

## -----------------------------------------------------------------------------
epane <- function(size){
  U1 <- runif(size, min = -1, max = 1)
  U2 <- runif(size, min = -1, max = 1)
  U3 <- runif(size, min = -1, max = 1)
  
  uu <- abs(U3) > abs(U1) & abs(U3) > abs(U2)
  x <- integer(size)
  x[uu] <- U2[uu]
  x[!uu] <- U3[!uu]
}

set.seed(123)
x <- epane(100000)

## ----echo=FALSE---------------------------------------------------------------
dens <- function(x) return(3 * (1-x^2) / 4)

p <- ggplot(data.frame(x), aes(x = x)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    colour = 1,
    fill = "gray"
    ) +
  stat_function(
    fun = dens,
    linewidth = 1
  ) + 
  ggtitle("样本频率与理论密度对比图") +
  theme(plot.title = element_text(hjust = 0.5))
p

## -----------------------------------------------------------------------------
set.seed(12345)
d <- 1; m <- 1e6
rho <- c(0.2, 0.5, 1)

pihat.var <- integer(0)
for(rr in rho){
  l <- d * rr
  pihat <- integer(0)
  for(i in 1:100){
    X <- runif(m,0,d/2)
    Y <- runif(m,0,pi/2)
    pihat <- c(pihat, 2*l/d/mean(l/2*sin(Y)>X))
  }
  pihat.var <- c(pihat.var, m * var(pihat))
}

pihat.var

## -----------------------------------------------------------------------------
pi^3 / 2 / rho - pi^2

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1e6
u <- runif(n)

theta.mc <- mean(exp(u))#简单MC

u <- u[1:(n/2)]#对偶变量法
theta.ant <- mean((exp(u) + exp(1-u)) / 2)

print(c(theta.mc, theta.ant))

## -----------------------------------------------------------------------------
set.seed(123)
n <- 1e6

var.mc <- integer(100)
var.ant <- integer(100)

for(i in 1:100){
  u <- runif(n)
  var.mc[i] <- var(exp(u) / n)#简单MC
  
  u <- u[1:(n/2)]#对偶变量法
  var.ant[i] <- var((exp(u) + exp(1-u)) / n)
}

mean(1 - var.ant / var.mc)

## -----------------------------------------------------------------------------
x <- seq(1, 10, 0.01)
g <- x^2 * exp(-x^2/2) / sqrt(2*pi)
f1 <- 2*exp(-(x-1)^2/2) / sqrt(2*pi)
f2 <- exp(-x+1)
plot(x, g, type = "l", ylim = c(0, 1))
lines(x, f1, lty = 2)
lines(x, f2, lty = 3)
legend(
  "topright",
  legend = c("g", "f1", "f2"), lty = 1:3
)

## -----------------------------------------------------------------------------
plot(x, g/f1, type = "l",
     ylim = c(0, 0.6), lty = 2, ylab = "")
lines(x, g/f2, lty = 3)
legend(
  "topright",
  legend = c("g/f1", "g/f2"), lty = 2:3
)

## -----------------------------------------------------------------------------
set.seed(123)

m <- 1e4
num.rep <- 1e3
mc <- numeric(num.rep)

for(i in 1:num.rep){
  x <- abs(rnorm(m, 0, 1)) + 1
  g <- x^2 * exp(-x^2/2) / sqrt(2*pi)
  f <- 2*exp(-(x-1)^2/2) / sqrt(2*pi)
  
  mc[i] <- mean(g/f)
}

c(mean(mc), var(mc))

## -----------------------------------------------------------------------------
g <- function(x) x^2 * exp(-x^2/2) / sqrt(2*pi)
integrate(g, lower = 1, upper = Inf)

## -----------------------------------------------------------------------------
set.seed(123)

M <- 1e4; k <- 5
xx <- seq(0, 1, 1/k)

theta.hat <- theta.var <- numeric(k)
for(i in 1:k){
  u <- runif(M/k, min = xx[i], max = xx[i+1])
  x <- -log(1 - (1 - exp(-1)) * u)
  
  f <- (k/(1-exp(-1))) * exp(-x)
  g <- exp(-x) / (1+x^2)
  theta.hat[i] <- mean(g/f)
  theta.var[i] <- var(g/f)
}

c(sum(theta.hat), mean(theta.var))

## -----------------------------------------------------------------------------
sqrt(mean(theta.var))

## -----------------------------------------------------------------------------
set.seed(123)

num.rep <- 1e4; n <- 20
t0 <- qt(c(0.05/2, 1 - 0.05/2), df = n-1)

tt <- matrix(NA, nrow = num.rep, ncol = 2)
for(i in 1:num.rep){
  x <- rchisq(n, df = 2)
  tt[i, ] <- mean(x) + t0 * sd(x) / sqrt(n)
}

mean(tt[,1] < 2 & tt[,2] > 2)

## ----eval=FALSE---------------------------------------------------------------
#  # 1: chi方 2: 均匀 3: 指数
#  samp.generate <- function(num.rep, size){
#    samp <- array(NA, dim = c(num.rep, size, 3))
#    for(i in 1:num.rep){
#      samp[i, , 1] <- rchisq(size, df = 1)
#      samp[i, , 2] <- runif(size, min = 0, max = 2)
#      samp[i, , 3] <- rexp(size, rate = 1)
#    }
#    saveRDS(samp, file = "samp.rds")
#  }

## ----eval=FALSE---------------------------------------------------------------
#  samp.analysis <- function(filename, alpha){
#    samp <- readRDS(filename)
#    num.rep <- dim(samp)[1]; size <- dim(samp)[2]
#  
#    t0 <- qt(c(alpha/2, 1 - alpha/2), df = size-1)
#    tt <- array(NA, dim = c(num.rep, 2, 3))
#  
#    for(i in 1:3) for(j in 1:num.rep){
#      x <- samp[j, , i]
#      tt[j, , i] <- mean(x) + t0 * sd(x) / sqrt(size)
#    }
#    saveRDS(tt, file = paste(alpha, "interval.rds", sep = " "))
#  }

## ----eval=FALSE---------------------------------------------------------------
#  samp.res <- function(filename){
#    tt <- readRDS(filename)
#    num.rep <- dim(tt)[1]
#    is.in <- matrix(NA, nrow = num.rep, ncol = 3)
#    for(i in 1:3)
#      for(j in 1:num.rep){
#        is.in[j, i] <- tt[j, 1, i] < 1 & tt[j ,2, i] > 1
#      }
#    res <- c(1 - mean(is.in[,1]), 1 - mean(is.in[,2]), 1 - mean(is.in[,3]))
#    saveRDS(res, file = "result.rds")
#  }

## -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)
# 每次调用函数之前会清理所有内存。调用函数的代码已省略。

## ----echo=FALSE---------------------------------------------------------------
samp.generate <- function(num.rep, size){
  samp <- array(NA, dim = c(num.rep, size, 3))
  for(i in 1:num.rep){
    samp[i, , 1] <- rchisq(size, df = 1)
    samp[i, , 2] <- runif(size, min = 0, max = 2)
    samp[i, , 3] <- rexp(size, rate = 1)
  }
  saveRDS(samp, file = "samp.rds")
}

## -----------------------------------------------------------------------------
samp.generate(num.rep = 1e4, size = 1000)


## -----------------------------------------------------------------------------
rm(list = ls())

## ----echo=FALSE---------------------------------------------------------------
samp.analysis <- function(filename, alpha){
  samp <- readRDS(filename)
  num.rep <- dim(samp)[1]; size <- dim(samp)[2]
  
  t0 <- qt(c(alpha/2, 1 - alpha/2), df = size-1)
  tt <- array(NA, dim = c(num.rep, 2, 3))
  
  for(i in 1:3) for(j in 1:num.rep){
    x <- samp[j, , i]
    tt[j, , i] <- mean(x) + t0 * sd(x) / sqrt(size)
  }
  saveRDS(tt, file = paste(alpha, "interval.rds", sep = " "))
}

## -----------------------------------------------------------------------------
samp.analysis("samp.rds", alpha = 0.05)

rm(list = ls())

## ----echo=FALSE---------------------------------------------------------------
samp.res <- function(filename){
  tt <- readRDS(filename)
  num.rep <- dim(tt)[1]
  is.in <- matrix(NA, nrow = num.rep, ncol = 3)
  for(i in 1:3) 
    for(j in 1:num.rep){
      is.in[j, i] <- tt[j, 1, i] < 1 & tt[j ,2, i] > 1
    }
  res <- c(1 - mean(is.in[,1]), 1 - mean(is.in[,2]), 1 - mean(is.in[,3]))
  saveRDS(res, file = "result.rds")
}

## -----------------------------------------------------------------------------
samp.res("0.05 interval.rds")

##每个函数相互独立，只需要读取前一函数生成的文件。

result <- readRDS("result.rds")
result

## ----echo=FALSE---------------------------------------------------------------
samp.analysis <- function(filename, alpha){
  samp <- readRDS(filename)
  num.rep <- dim(samp)[1]; size <- dim(samp)[2]
  
  t0 <- qt(c(alpha/2, 1 - alpha/2), df = size-1)
  tt <- array(NA, dim = c(num.rep, 2, 3))
  
  for(i in 1:3) for(j in 1:num.rep){
    x <- samp[j, , i]
    tt[j, , i] <- mean(x) + t0 * sd(x) / sqrt(size)
  }
  saveRDS(tt, file = paste(alpha, "interval.rds", sep = " "))
}


samp.analysis("samp.rds", alpha = 0.1)
samp.res("0.1 interval.rds")
result <- readRDS("result.rds")
result

## -----------------------------------------------------------------------------
#
# bon: bonferroni bh: Benjamini-Hochberg
#
set.seed(123)
m <- 1000; M <- 1000

bon.fwer <- bh.fwer <- numeric(M)
bon.fdr <- bh.fdr <- numeric(M)
bon.tpr <- bh.tpr <- numeric(M)

for(i in 1:M){
  p <- c(runif(0.95 * m), rbeta(0.05 * m, 0.1, 1))
  
  p.bon <- p.adjust(p, method = "bonferroni")
  ord <- order(p)
  p.bh <- p.adjust(p[ord], method = "fdr")
  p.bh[ord] <- p.bh
  
  bon.fwer[i] <- sum(head(p.bon, 0.95 * m) < 0.1) > 0
  bh.fwer[i] <- sum(head(p.bh, 0.95 * m) < 0.1) > 0
  
  bon.fdr[i] <- sum(head(p.bon, 0.95 * m) < 0.1) / sum(p.bon < 0.1)
  bh.fdr[i] <- sum(head(p.bh, 0.95 * m) < 0.1) / sum(p.bh < 0.1)
  
  bon.tpr[i] <- mean(tail(p.bon, 0.05 * m) < 0.1)
  bh.tpr[i] <- mean(tail(p.bh, 0.05 * m) < 0.1)
}

result <- matrix(NA, nrow = 2, ncol = 4)
result[1,] <- c("Bonferroni", round(c(mean(bon.fwer), mean(bon.fdr), mean(bon.tpr)), 5))
result[2,] <- c("B-H", round(c(mean(bh.fwer), mean(bh.fdr), mean(bh.tpr)), 5))
colnames(result) <- c("method", "FWER", "FDR", "TPR")


knitr::kable(result, align = "c")

## -----------------------------------------------------------------------------
rm(list = ls())
set.seed(123)

lam <- 2; nn <- c(5, 10, 20)
m <- 1000; B <- 1000

result <- matrix(NA, nrow = 3, ncol = 5)
colnames(result) <- c("size", "bias-bs", "se-bs", "bias-th", "se-th")
result[, 1] <- nn

for(i in 1:3){
  n <- nn[i]
  bias <- se <- numeric(m)
  
  for(mm in 1:m){
    samp <- rexp(n, rate = lam)
    
    lam.hat <- numeric(m)
    for(b in 1:B){
      samp.b <- sample(samp, n, replace = TRUE)
      lam.hat[b] <- 1 / mean(samp.b)
    }
    bias[mm] <- mean(lam.hat) - lam
    se[mm] <- sd(lam.hat)
  }
  result[i, 2:3] <- round(c(mean(bias), mean(se)), 5)
  result[i, 4:5] <- round(c(lam / (n-1), lam*n / (n-1) / sqrt(n-2)), 5)
}

#
#bs: bootstrap th: 理论值
#
knitr::kable(result, align = "c")

## -----------------------------------------------------------------------------
rm(list = ls())
law <- bootstrap::law
x <- law$LSAT; y <- law$GPA

B <- 1000; BB <- 50 #BB: 使用bootstrap计算se的重复次数
alpha <- 0.05

set.seed(123)
n <- length(x)
cor.obs <- cor(x, y)

t <- cor.bs <- numeric(B)
for(b in 1:B){
  samp <- sample(1:n, n, replace = TRUE)
  x.samp <- x[samp]; y.samp <- y[samp]
  cor.bs[b] <- cor(x.samp, y.samp)
  
  cor.rebs <- numeric(BB)
  for(bb in 1:BB){
    resamp <- sample(1:n, n, replace = TRUE)
    cor.rebs[bb] <- cor(x.samp[resamp], y.samp[resamp])
  }
  t[b] <- (cor.bs[b] - cor.obs) / sd(cor.rebs)
}

CI <- mean(cor.bs) - quantile(t, c(1-alpha/2, alpha/2)) * sd(cor.bs)
names(CI) <- rev(names(CI))

CI

## -----------------------------------------------------------------------------
rm(list = ls())
library(boot)
set.seed(123)

f <- function(x, i) mean(x[i,])

b <- boot(data = aircondit, statistic = f, R = 1000)
CI <- boot.ci(b, conf = 0.95, type=c("norm", "basic", "perc", "bca"))
# CI

result <- matrix(NA, nrow = 4, ncol = 3)
colnames(result) <- c("method", "lower", "upper")

result[, 1] <- c("norm", "basic", "perc", "bca")
result[1, 2:3] <- round(CI$normal[2:3], 3); result[2, 2:3] <- round(CI$basic[4:5], 3)
result[3, 2:3] <- round(CI$perc[4:5], 3); result[4, 2:3] <- round(CI$bca[4:5], 3)

knitr::kable(result, align = "c")

## -----------------------------------------------------------------------------
hist(b$t, prob = TRUE)
abline(v = b$t0)

## -----------------------------------------------------------------------------
rm(list = ls())
library(bootstrap)

n <- nrow(scor)

theta.est <- function(x){
  #估计theta
  v <- eigen(cov(x))$values
  return(max(v) / sum(v))
}

theta.hat <- theta.est(scor)
theta.jack <- sapply(1:n, FUN = function(i) theta.est(scor[-i, ]))

bias <- (n-1) * (mean(theta.jack) - theta.hat)
se <- sqrt(n-1) * sd(theta.jack)

result <- c(theta.hat, bias, se)
names(result) <- c("est", "bias", "se")
result

## -----------------------------------------------------------------------------
rm(list = ls())
library(DAAG); attach(ironslag)
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n * (n-1) / 2)

fff <- 1
for(k in 1:(n-1)) for(l in (k+1):n){
  y <- magnetic[-c(k, l)]
  x <- chemical[-c(k, l)]
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[c(k, l)]
  e1[fff] <- sum((magnetic[c(k, l)] - yhat1)^2)
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[c(k, l)] +
    J2$coef[3] * chemical[c(k, l)]^2
  e2[fff] <- sum((magnetic[c(k, l)] - yhat2)^2)
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[c(k, l)]
  yhat3 <- exp(logyhat3)
  e3[fff] <- sum((magnetic[c(k, l)] - yhat3)^2)
  
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[c(k, l)])
  yhat4 <- exp(logyhat4)
  e4[fff] <- sum((magnetic[c(k, l)] - yhat4)^2)
  
  fff <- fff + 1
}

c(mean(e1), mean(e2), mean(e3), mean(e4))

## -----------------------------------------------------------------------------
rm(list = ls())

cvm.test <- function(x, y){
  m <- length(x); n <- length(y)
  F1 <- ecdf(x); F2 <- ecdf(y)
  w <- sum((F1(x) - F2(x))^2) + sum((F1(y) - F2(y))^2)
  w <- m*n / (m+n)^2 * w
  return(w)
}

cvm.perm <- function(x, y, B){
  z <- c(x, y)
  perm <- numeric(B + 1)
  t0 <- cvm.test(x, y)
  
  m <- length(x); n <- length(y)
  
  for(i in 1:B){
    xy <- sample(z)
    x1 <- xy[1:m]; y1 <- xy[-(1:m)]
    perm[i] <- cvm.test(x1, y1)
  }
  perm[B + 1] <- t0
  
  p <- mean(perm >= t0)
  return(list(t0 = t0, perm = perm, p = p))
}

## -----------------------------------------------------------------------------
B <- 999
attach(chickwts)
x1 <- as.vector(weight[feed == "soybean"])
x2 <- as.vector(weight[feed == "linseed"])
x3 <- as.vector(weight[feed == "sunflower"])
detach(chickwts)

set.seed(123)

s1 <- cvm.perm(x1, x2, B)
s2 <- cvm.perm(x2, x3, B)

## -----------------------------------------------------------------------------
s1$p

## -----------------------------------------------------------------------------
s2$p

## -----------------------------------------------------------------------------
max.ext <- function(x, y){
  x <- x - mean(x); y <- y - mean(y)
  mx <- sum(x > max(y)) + sum(x < min(y))
  my <- sum(y > max(x)) + sum(y < min(x))
  return(max(c(mx, my)))
}

max.ext.perm <- function(x, y, B){
  m <- length(x); n <- length(y)
  s0 <- max.ext(x, y)
  z <- c(x, y)
  
  perm <- numeric(B + 1)
  
  for(i in 1:B){
    xy <- sample(z)
    x1 <- xy[1:m]; y1 <- xy[-(1:m)]
    perm[i] <- max.ext(x1, y1)
  }
  perm[B + 1] <- s0
  
  p <- mean(perm >= s0)
  
  return(list(s0 = s0, perm = perm, p = p))
}

## -----------------------------------------------------------------------------
set.seed(123)

x <- rnorm(50)
y <- rnorm(100)

res1 <- max.ext.perm(x, y, 999)
res1$p

## -----------------------------------------------------------------------------
set.seed(456)

x <- rnorm(50)
y <- rnorm(100, sd = 2)

res2 <- max.ext.perm(x, y, 999)
res2$p

## -----------------------------------------------------------------------------
rm(list = ls())

a.compute <- function(N, b, f0){
  X1 <- rpois(N, lambda = 1); X2 <- rexp(N); X3 <- rbinom(N, 1, 0.5)
  
  g <- function(a){
    tt <- exp(- a - b[1] * X1 - b[2] * X2 - b[3] * X3)
    return(mean(1 / (1 + tt)) - f0)
  }
  solution <- uniroot(g,c(-20,0))
  return(solution)
}

## -----------------------------------------------------------------------------
set.seed(123)
f0 <- c(0.1, 0.01, 0.001, 0.0001)
N <- 1e6; b <- c(0, 2, -1)

res <- lapply(f0, a.compute, N = N, b = b)

a <- sapply(res, function(u) u$root)
result <- rbind(f0 = f0, a = a)
knitr::kable(result, align = "c")

## -----------------------------------------------------------------------------
plot(log(f0), a)

## -----------------------------------------------------------------------------
rm(list = ls())

mh.Laplace <- function(N, x0, sigma){
  #N: 链长度| x0: 初值| sigma: 正态分布标准差
  x <- numeric(N + 1); x[1] <- x0
  u <- runif(N)
  num.acc <- 0
  
  for(i in 2:(N+1)){
    y <- rnorm(1, x[i-1], sigma)
    if(u[i-1] <= exp(abs(x[i-1]) - abs(y))){
      x[i] <- y
      num.acc <- num.acc + 1
    }else x[i] <- x[i-1]
  }
  
  return(list(chain = x, num.acc = num.acc))
}

## -----------------------------------------------------------------------------
set.seed(123)
N <- 5000; x0 <- rnorm(1)
sigma <- sqrt(c(0.1, 0.5, 1, 2, 4, 16))

res <- lapply(sigma, mh.Laplace, N = N, x0 = x0)
accept.rate <- sapply(1:6, function(i) res[[i]]$num.acc / N)

result <- rbind(var = sigma^2, accept.rate)
knitr::kable(result, align = "c")

## -----------------------------------------------------------------------------
par(mfrow = c(2, 3))
for(i in 1:6) plot(
    res[[i]]$chain, ylab = "", type = "l",
    xlab = paste("var = ", (sigma^2)[i])
  )

## -----------------------------------------------------------------------------
lam.update <- function(lam, u, v){
  uu <- exp(-lam * u) - exp(-lam * v)
  ll <- u * exp(-lam * u) - v * exp(-lam * v)
  
  ss <- 1 / lam + mean(ll / uu)
  
  return(1 / ss)
}

lam.em <- function(init, u, v, tol = 1e-5){
  est.old <- init
  est <- lam.update(init, u, v)
  
  while(abs(est.old - est) > tol){
    est.old <- est
    est <- lam.update(est, u, v)
  }
  return(est)
}

#(11,12), (8,9), (27,28), (13,14), (16,17), (0,1), (23,24), (10,11), (24,25), (2,3)
u <- c(11, 8, 27, 13, 16, 0, 23, 10, 24, 2)
v <- c(12, 9, 28, 14, 17, 1, 24, 11, 25, 3)

MLE1 <- lam.em(1, u, v)
MLE1

## -----------------------------------------------------------------------------
f <- function(lam, u, v){
  s <- u * exp(-lam * u) - v * exp(-lam *v)
  t <- exp(-lam * u) - exp(-lam *v)
  return(sum(s / t))
}

MLE2 <- uniroot(f, c(1e-3, 1), u = u, v = v, tol = 1e-5)
MLE2$root

## -----------------------------------------------------------------------------
solve.game <- function(A){
  min.A <- min(A)
  A <- A - min.A
  max.A <- max(A)
  A <- A / max(A)
  m <- nrow(A)
  n <- ncol(A)
  it <- n^3
  a <- c(rep(0, m), 1)
  A1 <- -cbind(t(A), rep(-1, n))
  b1 <- rep(0, n)
  A3 <- t(as.matrix(c(rep(1, m), 0)))
  b3 <- 1
  sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
      maxi=TRUE, n.iter=it)
  a <- c(rep(0, n), 1)
  A1 <- cbind(A, rep(-1, m))
  b1 <- rep(0, m)
  A3 <- t(as.matrix(c(rep(1, n), 0)))
  b3 <- 1
  sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
  maxi=FALSE, n.iter=it)
  soln <- list("A" = A * max.A + min.A,
  "x" = sx$soln[1:m],
  "y" = sy$soln[1:n],
  "v" = sx$soln[m+1] * max.A + min.A)
  soln
}

A <- matrix(c(
  0,-2,-2,3,0,0,4,0,0,
  2,0,0,0,-3,-3,4,0,0,
  2,0,0,3,0,0,0,-4,-4,
  -3,0,-3,0,4,0,0,5,0,
  0,3,0,-4,0,-4,0,5,0,
  0,3,0,0,4,0,-5,0,-5,
  -4,-4,0,0,0,5,0,0,6,
  0,0,4,-5,-5,0,0,0,6,
  0,0,4,0,0,5,-6,-6,0), 9, 9)

library(boot)
B <- A + 2
s <- solve.game(B)

## -----------------------------------------------------------------------------
round(cbind(s$x, s$y), 5)

## -----------------------------------------------------------------------------
s$v

## -----------------------------------------------------------------------------
a <- c(1, 2, 3)
dim(a)

## -----------------------------------------------------------------------------
x <- matrix(0, 2, 2)
c(is.matrix(x), is.array(x))

## -----------------------------------------------------------------------------
m <- matrix(0, 2, 2)
m <- m[FALSE, FALSE]
m

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
x <- as.data.frame(matrix(rnorm(9), 3, 3)); x
y <- apply(x, 2, scale01); y

## -----------------------------------------------------------------------------
x <- cbind(x, s = c("a", "b", "c")); x
y <- apply(x, 2, function(u) if(is.numeric(u)) scale01(u) else u); y

## -----------------------------------------------------------------------------
x <- as.data.frame(matrix(rnorm(9), 3, 3)); x
y <- vapply(x, sd, numeric(1)); y

## -----------------------------------------------------------------------------
x <- cbind(x, s = c("a", "b", "c")); x
y <- vapply(
  x[vapply(x, is.numeric, logical(1))],
  sd, numeric(1)
)
y

## -----------------------------------------------------------------------------
gibbs_r <- function(a, b, n, N){
  X <- matrix(0, N, 2)
  
  for(i in 2:N){
    X[i, 1] <- rbinom(1, n, X[i-1, 2])
    X[i, 2] <- rbeta(1, X[i, 1] + a, n - X[i, 1] + b)
  }
  
  colnames(X) <- c("x", "y")
  return(X)
}

## -----------------------------------------------------------------------------
library(Rcpp)

cppFunction('
NumericMatrix gibbs_rcpp(double a, double b, int n, int N){
  NumericMatrix X(N, 2);
  
  for(int i = 1; i < N; i++){
    X(i, 0) = R::rbinom(n, X(i-1, 1));
    X(i, 1) = R::rbeta(X(i, 0) + a, n - X(i, 0) + b);
  }
  
  return X;
}     
')

## -----------------------------------------------------------------------------
library(microbenchmark)
a <- 3; b <- 4; n <- 20; N <- 1e4

microbenchmark(gibbs_r(a, b, n, N), gibbs_rcpp(a, b, n, N))

