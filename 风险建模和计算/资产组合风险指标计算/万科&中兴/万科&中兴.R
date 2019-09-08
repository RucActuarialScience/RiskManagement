library(quantmod)
library(ggplot2)
library(mvnormtest)
library(QRM)
library(zoo)
library(mvtnorm) # for sampling from a multivariate normal or t distribution
library(rmgarch)
library(fBasics)
library(qrmtools)
library(forecast)
library(utils)
options(digits = 3)
#################GET THE DATA##########################################
set.seed(12)
setSymbolLookup(WK = list(name = '000002.sz', src = 'yahoo'))#万科
getSymbols("WK", from = '2011-01-01', to = '2015-12-31', warnings = FALSE)
#dailyReturn(WK)
WK = na.omit(WK)#去除NA数据
daily_return_WK = periodReturn(WK, period = 'daily') 
monthly_return_WK = periodReturn(WK, period = 'monthly') 
chartSeries(WK)

setSymbolLookup(ZXTX = list(name = '000063.sz', src = 'yahoo'))#中兴通讯
getSymbols("ZXTX", from = '2011-01-01', to = '2015-12-31')
ZXTX = na.omit(ZXTX)
daily_return_ZXTX = periodReturn(ZXTX, period = 'daily') 
monthly_return_ZXTX = periodReturn(ZXTX, period = 'monthly') 
chartSeries(ZXTX)


chartSeries(WK)
chartSeries(ZXTX)
plot(daily_return_ZXTX)
plot(monthly_return_ZXTX)
plot(daily_return_WK)
plot(monthly_return_WK)


#####################Tests of multivariate normality####################

dailyret = cbind(daily_return_WK, daily_return_ZXTX)
dailyret = na.omit(dailyret)
colnames(dailyret) = c("WK", "ZXTX")

monthlyret = cbind(monthly_return_WK, monthly_return_ZXTX)
monthlyret = na.omit(monthlyret)
colnames(monthlyret) = c("WK", "ZXTX")
plot.zoo(dailyret)
plot.zoo(monthlyret)

########daily returns
apply(dailyret, 2, function(x) shapiro.test(x)$p.value)
jointnormalTest(as.matrix(dailyret), plot = FALSE)
MardiaTest(as.matrix(dailyret))

########monthly returns
apply(monthlyret, 2, function(x) shapiro.test(x)$p.value)
jointnormalTest(as.matrix(monthlyret), plot = FALSE)
MardiaTest(as.matrix(monthlyret))

## Visual tests of normality

## Daily returns
pairs(as.matrix(dailyret), gap = 0, pch = ".") # visual assessment
par(mfrow = c(1, 2))
qq_plot(dailyret[, 1], FUN = qnorm, method = "empirical") # first margin only = > already not normal
qq_plot(dailyret[, 2], FUN = qnorm, method = "empirical")
D2.d = mahalanobis(dailyret, center = colMeans(dailyret), cov = cov(dailyret)) # squared Mahalanobis distances
par(mfrow = c(1, 1))
qq_plot(D2.d , FUN = function(p) qchisq(p, df = 2)) # = > departure clearly visible


## monthly returns
pairs(as.matrix(monthlyret), gap = 0, pch = ".") # visual assessment
par(mfrow = c(1, 2))
qq_plot(monthlyret[, 1], FUN = qnorm, method = "empirical") 
qq_plot(monthlyret[, 2], FUN = qnorm, method = "empirical") 
D2.m = mahalanobis(monthlyret, center = colMeans(monthlyret), cov = cov(monthlyret)) # squared Mahalanobis distances
par(mfrow = c(1, 1))
qq_plot(D2.m, FUN = function(p) qchisq(p, df = 2)) # = > departure clearly visible


###################### Fit various multivariate models###########################

## Fitting a multivariate normal distribution to X and simulating from it
mu = colMeans(dailyret) # estimated location vector
Sigma = cov(dailyret) # estimated scale matrix
stopifnot(all.equal(Sigma, var(dailyret)))
P = cor(dailyret) # estimated correlation matrix
stopifnot(all.equal(P, cov2cor(Sigma)))
n = nrow(dailyret) # sample size
set.seed(866)
X.norm = rmvnorm(n, mean = mu, sigma = Sigma) # N(mu, Sigma) samples

## Fitting a multivariate t distribution to X
fit = fit.mst(dailyret, method = "BFGS") # fit a multivariate t distribution
X.t = rmvt(n, sigma = as.matrix(fit$Sigma), df = fit$df, delta = fit$mu) # t_nu samples

## Plot (sample from fitted t (red), original sample (black), sample from fitted normal (blue))

pairs(as.matrix(dailyret), gap = 0, pch = ".", col = "black")
pairs(X.norm, gap = 0, pch = ".", col = "blue")
pairs(X.t, gap = 0, pch = ".", col = "red")

#par(mfrow = c(1 , 1))
#plot(as.matrix(dailyret)[, 1:2] , col = "black")
#plot(X.norm[, 1:2] , col = "blue")
#plot(X.t[, 1:2] , col = "red")
ggplot(data = dailyret, mapping = aes(WK, ZXTX))+geom_point()
ggplot(data = as.data.frame(X.norm), mapping = aes(WK, ZXTX))+geom_point()
colnames(X.t) = c("WK", "ZXTX")
ggplot(data = as.data.frame(X.t), mapping = aes(WK, ZXTX))+geom_point()

###########calculate VAR ES#####
alpha = 0.95
mu_new1 = 0.5*mu[1]+0.5*mu[2]
sigma_new1 = sqrt(0.25*Sigma[1, 1]^2+0.25*Sigma[2, 2]^2)
X.norm_new1 = rnorm(n, mu_new1, sigma_new1)
Var1 = quantile(X.norm_new1, probs = alpha)
ES1 = sum(X.norm_new1*(X.norm_new1>Var1))/(sum(X.norm_new1>Var1))

## Specify univariate AR(1)-GARCH(1, 1) for both marginal processes#########

acf(dailyret)

uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
     mean.model = list(armaOrder = c(1, 0), include.mean = TRUE), 
     distribution.model = "std")

## Check the univariate specification for the two component series
fit.marg1 <- ugarchfit(spec = uspec, data = dailyret[, 1])
fit.marg2 <- ugarchfit(spec = uspec, data = dailyret[, 2])

## Combine univariate specs to obtain spec for marginal models
marginspec <- multispec(replicate(2, uspec))

#############Create spec for DCC
mspec <- dccspec(marginspec, dccOrder = c(1, 1), model = "DCC", distribution = "mvt")

mod <- dccfit(mspec, dailyret)
mod
## Check marginal coefficients are same in joint model
coef(mod)
coef(fit.marg1)
coef(fit.marg2)
## Some pictures of fit

plot(mod, which = 2)
plot(mod, which = 3)
plot(mod, which = 4)
plot(mod, which = 5)

############VAR SE#####################
alpha = 0.95
ugarchforecast(fit.marg1, n.ahead = 1)
ugarchforecast(fit.marg2, n.ahead = 1)
mu_new2 = -0.5*(8.857e-05)-0.5*(9.664e-06)
sigma_new2 = sqrt(0.25*0.03483^2+0.25*0.02529^2)
X.norm_new2 = rnorm(n, mu_new2, sigma_new2)
Var2 = quantile(X.norm_new2, probs = alpha)
ES2 = sum(X.norm_new2*(X.norm_new2>Var2))/(sum(X.norm_new2>Var2))

####################将数据导出为EXCEL#########################
write.csv(WK, "E:/RUC/lessons/量化风险管理/第三次作业/数据/WK.csv")
write.csv(daily_return_WK, "E:/RUC/lessons/量化风险管理/第三次作业/数据/daily_return_WK.csv")
write.csv(monthly_return_WK, "E:/RUC/lessons/量化风险管理/第三次作业/数据/monthly_return_WK.csv")
write.csv(ZXTX, "E:/RUC/lessons/量化风险管理/第三次作业/数据/ZXTX.csv")
write.csv(daily_return_ZXTX, "E:/RUC/lessons/量化风险管理/第三次作业/数据/daily_return_ZXTX.csv")
write.csv(monthly_return_ZXTX, "E:/RUC/lessons/量化风险管理/第三次作业/数据/monthly_return_ZXTX.csv")