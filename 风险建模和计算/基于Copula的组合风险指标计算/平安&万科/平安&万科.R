library(MASS)
library(quantmod)
library(copula)
library(qrmtools)

#####1 Visualizing Sklar's Theorem
###1.1 绘制二元正态和高斯Copula样本散点图
set.seed(100)
sigma <- matrix(c(1, 0.3, 0.3, 1), ncol = 2)
mu <- c(0, 0)
sample.norm <- mvrnorm(n = 1000, mu, sigma)
par(mfrow = c(1,1), cex.main = 1)
plot(sample.norm[, 1], sample.norm[, 2], xlab = 'X', ylab = 'Y', 
     main = '1000 realizations of (X,Y) for a joint normal distribution with tho  = 0.3')
sample.copula <- pnorm(sample.norm)

plot(sample.copula, xlim = c(0, 1), ylim = c(0, 1),
     xlab = 'U = F(X)', ylab = 'V = F(Y)', main = '(U,V) = (F(X),F(Y)) for (X,Y) a joint normal distribution with tho  = 0.3')

###1.2 绘制相依结构为高斯copula，边缘分布为Exp(4)的样本散点图
lambda <- 4
sample.exp <- - log(1 - sample.copula) / lambda
plot(sample.exp, xlab = expression(G^-1~(U)), ylab = expression(H^-1~(V)),
     main = expression('('~G^-1~(U)~','~H^-1~(V)~') for dfs G, H exponential with mean 0.25'), cex.lab=0.7)

#####2 研究第三次作业日收益率#####
# 下载数据
setSymbolLookup(PA = list(name = '000001.sz', src = 'yahoo'))
setSymbolLookup(WK = list(name = '000002.sz', src = 'yahoo'))
getSymbols('PA', from = '2014-01-01', to = '2019-01-01')
# 只留下收盘数据
PA <- PA[, 4]
WK <- WK[, 4]

# 计算日对数收益率
PA_Daily <- Delt(PA, type = 'log')[-1]
WK_Daily <- Delt(WK, type = 'log')[-1]
R.Daily <- merge.xts(PA_Daily, WK_Daily)
names(R.Daily) <- c("PA", "WK")
R.Daily <- as.matrix(R.Daily)
plot(R.Daily, xlab=expression(X[1]), ylab=expression(X[2]), xlim = c(-0.15, 0.1), ylim = c(-0.15, 0.1), main = '平安银行和万科日收益率')

###2.1 计算相关系数
compute.kendall <- function(mat){
  n <- nrow(mat)
  kendall <- 0
  for(i in 2:n){
    diff <- sign(mat[,] - mat[i, ])
    kendall = kendall + sum(diff[1:i-1,1] * diff[1:i-1,2])
  }
  return(kendall / (n * (n-1) / 2))
}
compute.spearman <- function(mat){
  n <- nrow(mat)
  df <- pobs((mat)) # 计算经验分布函数值
  return(cor(df)[1,2])
}
(pearson <- cov(R.Daily)[1,2] / (sd(R.Daily[,1]) * sd(R.Daily[,2])))
(kendall <- compute.kendall(R.Daily))
(spearman <- compute.spearman(R.Daily))

###2.2 绘制上尾相依图和下尾相依图
conditional.prob <-function(mat,  quantiles, lower.tail){
  prob <- c()
  for(i in 1:nrow(quantiles)){
    if(lower.tail){
      mask <- mat[,1] <= quantiles[i, 1]
      prob <- c(prob, mean(mat[mask, 2] <= quantiles[i, 2]))
    }
    else{
      mask <- mat[,1] > quantiles[i,1]
      prob <- c(prob, mean(mat[mask, 2] > quantiles[i, 2])) 
      }
  }
  return(prob)
}
# 上尾相依图
u <- seq(0.9, 0.999, length.out = 100)
quantiles <- cbind(quantile(R.Daily[,1], u), quantile(R.Daily[,2], u))  # 计算相应的分位数
upper <- conditional.prob(R.Daily, quantiles, FALSE)
plot(u, upper, ylab = '概率', main = '上尾相依图')

# 下尾相依图
u <- seq(0.001, 0.1, length.out = 100)
quantiles <- cbind(quantile(R.Daily[,1], u), quantile(R.Daily[,2], u))  # 计算相应的分位数
lower <- conditional.prob(R.Daily, quantiles, TRUE)
plot(u, lower, ylab = '概率', main = '下尾相依图')

#####3 对日收益率数据进行模型拟合#####

###3.1 分别用高斯分布和t分布拟合边缘分布
# 拟合正态分布
(fit1.norm <- fitdistr(R.Daily[,1], 'normal'))
(fit2.norm <- fitdistr(R.Daily[,2], 'normal'))

# 拟合t分布
mydt <- function(x, m, s, df) dt((x-m)/s, df)/s
(fit1.t <- fitdistr(R.Daily[,1], mydt, list(m=0, s=1, df=4)))
(fit2.t <- fitdistr(R.Daily[,2], mydt, list(m=0, s=1, df=4)))

###3.2 分别拟合Gumbelcopula和t copula
U.empirical <- pobs(R.Daily)

U.norm <- matrix(0, ncol=2, nrow=1219)
U.norm[,1] <- pnorm(R.Daily[,1], mean=fit1.norm$estimate['mean'], sd=fit1.norm$estimate['sd'])
U.norm[,2] <- pnorm(R.Daily[,2], mean=fit2.norm$estimate['mean'], sd=fit2.norm$estimate['sd'])

U.t <- matrix(0, ncol=2, nrow=1219)
U.t[,1] <- pt((R.Daily[,1] - fit1.t$estimate['m']) / fit1.t$estimate['s'], df=fit1.t$estimate['df'])
U.t[,2] <- pt((R.Daily[,2] - fit2.t$estimate['m']) / fit2.t$estimate['s'], df=fit2.t$estimate['df'])

plot(U.empirical, xlab = expression(U[1]), ylab = expression(U[2]), xlim = c(0, 1), ylim = c(0, 1), main = '基于经验分布的copula样本')
plot(U.norm, xlab = expression(U[1]), ylab = expression(U[2]), xlim = c(0, 1), ylim = c(0, 1), main = '基于边缘正态分布的copula样本')
plot(U.t, xlab = expression(U[1]), ylab = expression(U[2]), xlim = c(0, 1), ylim = c(0, 1), main = '基于边缘t分布copula样本')

(fit.gumbel.emp <- fitCopula(gumbelCopula(), data=U.empirical))
(fit.t.emp <- fitCopula(tCopula(), data=U.empirical))

(fit.gumbel.norm <- fitCopula(gumbelCopula(), data=U.norm))
(fit.t.norm <- fitCopula(tCopula(), data=U.norm))

(fit.gumbel.t <- fitCopula(gumbelCopula(), data=U.t))
(fit.t.t <- fitCopula(tCopula(), data=U.t))

#####4 考虑四种情形下的VaR指标和ES指标#####
set.seed(100)
U.gumbel <- rCopula(1000, copula=fit.gumbel.emp@copula)
U.t <- rCopula(1000, copula=fit.t.emp@copula)
plot(U.gumbel, xlab=expression(U[1]), ylab=expression(U[2]), main='Gumbel copula sample')
plot(U.t, xlab=expression(U[1]), ylab=expression(U[2]), main='t copula sample')

# Meta Gumbel with Gauss margins
X1 <- matrix(0, ncol=2, nrow=1000)
X1[,1] <- qnorm(U.gumbel[,1], mean=fit1.norm$estimate['mean'], sd=fit1.norm$estimate['sd'])
X1[,2] <- qnorm(U.gumbel[,2], mean=fit2.norm$estimate['mean'], sd=fit2.norm$estimate['sd'])
plot(X1, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1), xlab=expression(X[1]), ylab=expression(X[2]), main='Meta Gumbel with Gauss margins')
Y1 <- rowMeans(X1)
VaR_np(Y1, 0.95)
ES_np(Y1, 0.95)

# Meta Gumbel with t margins
X2 <- matrix(0, ncol=2, nrow=1000)
X2[,1] <- qt(U.gumbel[,1], df=fit1.t$estimate['df']) * fit1.t$estimate['s'] + fit1.t$estimate['m']
X2[,2] <- qt(U.gumbel[,2], df=fit2.t$estimate['df']) * fit1.t$estimate['s'] + fit2.t$estimate['m']
plot(X2, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1), xlab=expression(X[1]), ylab=expression(X[2]), main='Meta Gumbel with t margins')
Y2 <- rowMeans(X2)
VaR_np(Y2, 0.95)
ES_np(Y2, 0.95)

# Meta t with Gauss margins
X3 <- matrix(0, ncol=2, nrow=1000)
X3[,1] <- qnorm(U.t[,1], mean=fit1.norm$estimate['mean'], sd=fit1.norm$estimate['sd'])
X3[,2] <- qnorm(U.t[,2], mean=fit2.norm$estimate['mean'], sd=fit2.norm$estimate['sd'])
plot(X3, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1), xlab=expression(X[1]), ylab=expression(X[2]), main='Meta t with Gauss margins')
Y3 <- rowMeans(X3)
VaR_np(Y3, 0.95)
ES_np(Y3, 0.95)

# Meta t with t margins
X4 <- matrix(0, ncol=2, nrow=1000)
X4[,1] <- qt(U.t[,1], df=fit1.t$estimate['df']) * fit1.t$estimate['s'] + fit1.t$estimate['m']
X4[,2] <- qt(U.t[,2], df=fit2.t$estimate['df']) * fit1.t$estimate['s'] + fit2.t$estimate['m']
plot(X4, xlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1), xlab=expression(X[1]), ylab=expression(X[2]), main='Meta t with t margins')
Y4 <- rowMeans(X4)
VaR_np(Y4, 0.95)
ES_np(Y4, 0.95)
