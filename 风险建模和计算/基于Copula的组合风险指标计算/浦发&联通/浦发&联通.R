library(mvtnorm)
library(MASS)
library(copula)
library(mvtnorm)
library(QRM)
## Sample from a norm copula

####### 1 ####################
n <- 1000 # sample size
d <- 2 # dimension
rho <- 0.3 # off-diagonal entry in the correlation matrix P
P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
set.seed(271)
X <- mvrnorm(n, c(rep(0, d)), P)
nc <- normalCopula(rho)
U.nc <- rCopula(n, copula = nc)
#X.nc <- qnorm(U.nc)
plot(X, xlab = quote(X[1]), ylab = quote(X[2])) #
plot(U.nc, xlab = quote(U[1]), ylab = quote(U[2]))#G-C
Y <- qexp(U.nc, rate = 4) # transform U.nc (Gauss copula) to Exp(4) margins
plot(Y, xlab = quote(Y[1]), ylab = quote(Y[2]))

#### 2 ############
library(quantmod)
#
setSymbolLookup(PF = list(name = '600000.ss', src = 'yahoo'))
getSymbols("PF")
#
setSymbolLookup(ZGLT = list(name = '600050.ss', src = 'yahoo'))
getSymbols("ZGLT")
G1 <- PF['2014-04-18/2019-05-29', 4]
G2 <- ZGLT['2014-04-18/2019-05-29', 4]
G <- cbind(G1, G2)
colnames(G) <- c("PF", "ZGLT")
## Daily log-returns
G.d <- returns(G)
G.d <- G.d[-1, ]
cor(G.d, method = "pearson")
cor(G.d, method = "kendall")
cor(G.d, method = "spearman")

#####
u = c(seq(0.9, 1, by = 0.001))
u = u[-length(u)]
p1 = NULL
for (i in 1:length(u)) {
  new = c(which(G.d[, 1] > quantile(ecdf(G.d[, 1]), u[i])))
  G.dd = G.d[new, ]
  p1[i] = length(which(G.dd[, 2] > quantile(ecdf(G.d[, 2]), u[i])))/length(G.dd[, 2])
}
#length(which(G.d[,1]>quantile(ecdf(G.d[,1]),u)))/length(G.d[,1])

### another way 
hh1 = c(which(G.d[, 1] > quantile(ecdf(G.d[, 1]), u)))
hh2 = c(which(G.d[, 2] > quantile(ecdf(G.d[, 2]), u)))
num1 = c(rep(0, length(G.d[, 1])))
for (i in 1:length(G.d[, 1])) {
  num1[hh1[i]]=1
}
num2 = c(rep(0, length(G.d[, 2])))
for (i in 1:length(G.d[, 2])) {
  num2[hh2[i]]=1
}
num = num1*num2
length(which(num > 0))/length(which(num1 > 0)) #right!
### end

plot(u, p1, type = "l")

####
u2 = c(seq(0, 0.1, by = 0.001))
u2 = u2[-1]
p2 = NULL
for (i in 1:length(u2)) {
  new = c(which(G.d[, 1] < quantile(ecdf(G.d[, 1]), u2[i])))
  G.dd = G.d[new, ]
  p2[i] = length(which(G.dd[, 2] < quantile(ecdf(G.d[, 2]), u2[i])))/length(G.dd[, 2])
}
plot(u2, p2, type = "l")

###### 3 ############
# 3.1
## Fitting a normal distribution to G.d and simulating from it
mu1 <- mean(G.d[, 1]) # estimated location vector
Sigma1 <- var(G.d[, 1]) # estimated scale matrix
mu2 <- mean(G.d[, 2]) # estimated location vector
Sigma2 <- var(G.d[, 2]) # estimated scale matrix

## Fitting a multivariate t distribution to X
fit1 <- fit.mst(G.d[, 1], method = "BFGS") # fit a multivariate t distribution
fit1$df
fit1$mu
fit1$gamma
fit1$Sigma
fit2 <- fit.mst(G.d[, 2], method = "BFGS") # fit a multivariate t distribution
fit2$df
fit2$mu
fit2$gamma
fit2$Sigma
# 3.2
U <- as.matrix(pobs(G.d))
fit.t <- fitCopula(tCopula(), data = U) # df of freedom are estimated, too
fit.G <- fitCopula(gumbelCopula(), data = U) #alpha=1.33

###### 4. ###########################
n=dim(G.d)[1]
th.g <- 1.334 # Gumbel copula parameter
gc <- gumbelCopula(th.g) # Gumbel copula
th.t <- 0.4086
tc <- tCopula(th.t) 

## Generate copula data
U.gc <- rCopula(n, copula = gc)
U.tc <- rCopula(n, copula = tc)

## Map to N(0,1) margins (meta-copula data)
X.gc <- qnorm(U.gc, sd = c(0.017, 0.027))
X.tc <- qnorm(U.tc, sd = c(0.017, 0.027))

#Map to t margins (meta-copula data)
T.gc <- qt(U.gc, c(1.961, 1.973))
T.tc <- qt(U.tc, c(1.961, 1.973))

#case 1
r1 = 0.5*X.gc[, 1] + 0.5*X.gc[, 2]
vaR1 = -quantile(r1, probs = 0.05)
cvaR1 = -mean(r1[which(r1 <= -vaR1)])
#case 2
r2 = 0.5*T.gc[, 1] + 0.5*T.gc[, 2]
vaR2 = -quantile(r2, probs = 0.05)
cvaR2 = -mean(r2[which(r2 <= -vaR2)])
#case 3
r3 = 0.5*X.tc[, 1] + 0.5*X.tc[, 2]
vaR3 = -quantile(r3, probs = 0.05)
cvaR3 = -mean(r3[which(r3 <= -vaR3)])
#case 4
r4 = 0.5*T.tc[, 1] + 0.5*T.tc[, 2]
vaR4 = -quantile(r4, probs = 0.05)
cvaR4 = -mean(r4[which(r4 <= -vaR4)])
