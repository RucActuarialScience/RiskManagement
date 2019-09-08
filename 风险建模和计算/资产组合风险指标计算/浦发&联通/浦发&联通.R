library(quantmod)
library(xts)
library(mvtnorm)
library(QRM)
library(qrmtools)
library(qrmdata)
library(ghyp)
library(rmgarch)
library(graphics)
library(lgarch)
#####1. input data #################
#浦发银行
setSymbolLookup(PF = list(name = '600000.ss', src = 'yahoo'))
getSymbols("PF")
#中国联通
setSymbolLookup(ZGLT = list(name = '600050.ss', src = 'yahoo'))
getSymbols("ZGLT")
G1 <- PF['2014-04-19/2019-05-04', 4]
G2 <- ZGLT['2014-04-19/2019-05-04', 4]
G <- cbind(G1, G2)
colnames(G) <- c("PF", "ZGLT")
## Daily log-returns
G.d <- returns(G)
plot.zoo(G.d)

## monthly log-returns
G.m <- apply.monthly(G.d, FUN = colSums)
plot.zoo(G.m)

###########  2. #############
# 2.1 p238 table 
MardiaTest(G.d) #daily b2 p value  k2 p value 
MardiaTest(G.m) #monthly b2 p value  k2 p value 

# 2.2 QQ plot
qq_plot(G.d[, 1], FUN = qnorm, method = "empirical") # first margin
qq_plot(G.d[, 2], FUN = qnorm, method = "empirical") # second margin

qq_plot(G.m[, 1], FUN = qnorm, method = "empirical") # first margin
qq_plot(G.m[, 2], FUN = qnorm, method = "empirical") # second margin

apply(G.d, 2, function(x) shapiro.test(x)$p.value) #检验边际分布正态性
apply(G.m, 2, function(x) shapiro.test(x)$p.value)

jointnormalTest(G.d) #daily QQ plot  squared Mahalanobis distances
jointnormalTest(G.m) #monthly qq plot squared Mahalanobis distances

#increase sample size only monthly 

GG1 <- PF['2008-04-19/',4]
GG2 <- ZGLT['2008-04-19/',4]
GG <- cbind(GG1,GG2)
colnames(GG) <- c("PF", "ZGLT")
## Daily log-returns
GG.d <- returns(GG)
GG.m <- apply.monthly(GG.d, FUN = colSums)

qq_plot(GG.m[,1], FUN = qnorm, method = "empirical") # first margin
qq_plot(GG.m[,2], FUN = qnorm, method = "empirical") # second margin

apply(GG.m, 2, function(x) shapiro.test(x)$p.value)
GG.m = ifelse(is.na(GG.m), 0,GG.m)
jointnormalTest(GG.m) #Reject the null hypothesis

########### 3. fit ##################

#3.1

## Fitting a multivariate normal distribution to G.d and simulating from it
mu <- colMeans(G.d) # estimated location vector
Sigma <- cov(G.d) # estimated scale matrix

stopifnot(all.equal(Sigma, var(G.d)))
P <- cor(G.d) # estimated correlation matrix
stopifnot(all.equal(P, cov2cor(Sigma)))
n <- nrow(G.d) # sample size
set.seed(271)
G.d.norm <- rmvnorm(n, mean = mu, sigma = Sigma) # N(mu, Sigma) samples

## Fitting a multivariate t distribution to X
fit <- fit.mst(G.d, method = "BFGS") # fit a multivariate t distribution
G.d.t <- rmvt(n, sigma = as.matrix(fit$Sigma), df = fit$df, delta = fit$mu) # t_nu samples

#3.2

## Plot (sample from fitted t (red), original sample (black), sample from fitted normal (blue))
dat1 <- rbind(t = G.d.t, original = as.matrix(G.d))
dat2 <- rbind(norm = G.d.norm, original = as.matrix(G.d))
cols <- rep(c("maroon3", "black"), each = n)

pairs(dat1, gap = 0, pch = ".", col = cols)
pairs(dat2, gap = 0, pch = ".", col = cols)
## Pick out one pair (to better see that multivariate t fits better)
plot(dat1, col = cols) # => the multivariate normal generates too few extreme losses!
legend("bottomright", bty = "n", pch = rep(1, 2), col = c("black", "maroon3"),
       legend = c("-Log-returns", expression("fitted" ~ italic(t)[nu](mu, Sigma))))

plot(dat2, col = cols) # => the multivariate normal generates too few extreme losses!
legend("bottomright", bty = "n", pch = rep(1, 2), col = c("black", "maroon3"),
       legend = c("-Log-returns", expression("fitted" ~ N(mu, Sigma))))

## compute VaR and ES
mu
Sigma

alpha <- qnorm(0.95, mean=0, sd=1)
VaR_0.95 = alpha*sqrt(0.25*(Sigma[1, 1] + Sigma[2, 2] + Sigma[1, 2] + Sigma[2,1])) - 0.25*(mu[1] + mu[2])
fa = exp(-alpha^2/2)/sqrt(2*pi)
ES_0.95 = sqrt(0.25*(Sigma[1, 1]+Sigma[2, 2]+Sigma[1, 2]+Sigma[2, 1]))*fa/0.05 + 0.25*(mu[1] + mu[2])

#### 4. Garch model #####################

#4.1

acf(G.d) #ZGLT和PF有延迟相关性，存在滞后关系，说明PF领先于ZGLT市场

#4.2 CCC-GARCH(1,1)
mod.ccc <- mlgarch(G.d)

print(mod.ccc)

##extract ccc-log-garch coefficients:
coef(mod.ccc)
##extract Gaussian log-likelihood (zeros excluded) of the ccc-log-garch model:
logLik(mod.ccc)
##extract Gaussian log-likelihood (zeros excluded) of the varma representation:
logLik(mod.ccc, varma = TRUE)
##extract variance-covariance matrix:
vcov(mod.ccc)

# DCC-GARCH(1,1)
uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
                    distribution.model = "norm")

## Check the univariate specification for the two component series
fit.marg1 <- ugarchfit(spec = uspec, data = G.d[,1])
fit.marg2 <- ugarchfit(spec = uspec, data =  G.d[,2])

## Combine univariate specs to obtain spec for marginal models
marginspec <- multispec(replicate(2, uspec))


## Create spec for DCC
mspec <- dccspec(marginspec, dccOrder = c(1,1), model = "DCC", distribution = "mvt")

mod <- dccfit(mspec,G.d)
print(mod)

## Check marginal coefficients are same in joint model
cof.dcc = coef(mod)
cof1 = coef(fit.marg1)
cof2 = coef(fit.marg2)

## Some pictures of fit
plot(mod, which = 3)
plot(mod, which = 4)
plot(mod, which = 5)


## VaR and ES
forc1 <- ugarchforecast(fit.marg1, n.ahead = 1)
forc2 <- ugarchforecast(fit.marg2, n.ahead = 1)

forc <- dccforecast(mod, n.ahead = 1)

alpha <- qnorm(0.95,mean=0,sd=1)
VaR.dcc_0.95 = alpha*sqrt(0.25*(0.01558^2 + 0.02643^2 + 2*0.0004117794*0.5811)) - 0.25*(-4.829e-05 + 0.0006352)
fa = exp(-alpha^2/2)/sqrt(2*pi)
ES.dcc_0.95 = sqrt(0.25*(0.01558^2 + 0.02643^2 + 2*0.0004117794*0.5811))*fa/0.05 + 0.25*(-4.829e-05 + 0.0006352)

