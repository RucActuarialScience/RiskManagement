library(zoo)
library(xts)
library(TTR)
library(quantmod)
# 下载股票数据
set.seed(271)
getSymbols("AAPL", from = "2003-01-01", to = Sys.Date(), src = "yahoo")
# 计算对数日收益率
library(PerformanceAnalytics)
rate <- periodReturn(AAPL$AAPL.Adjusted, period = "daily", type = "log")



## Q1
# 损失序列
loss.AAPL <- -rate
# 损失序列的时序图
chartSeries(
  loss.AAPL, type = "1", TA = NULL,
  subset = "2000/2019",
  name = "AAPL Loss",
  theme = "white", major.ticks = "years", minor.ticks = FALSE
)
# 损失序列的正态QQ图
x <- coredata(loss.AAPL)[,1]
qqnorm(x, main = "AAPL Loss")
qqline(x, col = "red")
# 利用Jarque-Bera进行正态性检验
library(fBasics) 
basicStats(x)
normalTest(x, method = "jb")
# 进行平稳性检验，确定ARMA模型
par(mfrow = c(1,2))
acf(loss.AAPL)
pacf(loss.AAPL)
library(TSA)
eacf(loss.AAPL)
resm1 <- arima(loss.AAPL, order = c(1,0,1))
resm1
resm2 <- arima(loss.AAPL, order = c(2,0,2))
resm2
resm3 <- arima(loss.AAPL, order = c(3,0,3))
resm3
Box.test(resm2$residuals, type = 'Ljung-Box')

# arma-garch模型拟合
acf(abs(loss.AAPL))
acf(loss.AAPL^2)
Box.test(abs(loss.AAPL), type = 'Ljung-Box')
Box.test(loss.AAPL^2, type = 'Ljung-Box')
library(rugarch)
uspec.N <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(2,2), include.mean = TRUE), distribution.model = "norm")
fit.N <- ugarchfit(spec = uspec.N, data = loss.AAPL)
layout(matrix(1:4, ncol = 2, byrow = TRUE))
plot(fit.N, which = 6) # ACF of |X_t|
plot(fit.N, which = 9) # Q-Q plot of Z_t against a normal
plot(fit.N, which = 10) # ACF of Z_t
plot(fit.N, which = 11) # ACF of Z_t^2
Box.test(fit.N@fit$residuals, type = "Ljung-Box")
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), mean.model = list(armaOrder = c(2,2), include.mean = TRUE), distribution.model = "std")
fit.t <- ugarchfit(spec = uspec.t, data = loss.AAPL)
plot(fit.t, which = 6)
plot(fit.t, which = 9)
plot(fit.t, which = 10)
plot(fit.t, which = 11)
Box.test(fit.t@fit$residuals, type = "Ljung-Box")
layout(1)
LL.N <- fit.N@fit$LLH
LL.t <- fit.t@fit$LLH
AIC.N <- 2*length(fit.N@fit$coef) - 2*LL.N
AIC.t <- 2*length(fit.t@fit$coef) - 2*LL.t
forc <- ugarchforecast(fit.t, n.ahead = 1)
forc

# 新息序列即为残差
zt <- fit.t@fit$residuals
zt <- merge(loss.AAPL, zt)[,2]
plot(zt)
head(zt)



## Q2
# 2.1
endpts <- endpoints(zt, "days")
endpts1 <- endpts[seq(1, length(endpts), by = 60)]
M.60 <- period.apply(zt, INDEX = endpts1, FUN = max)
endpts2 <- endpts[seq(1, length(endpts), by = 120)]
M.120 <- period.apply(zt, INDEX = endpts2, FUN = max)
endpts3 <- endpts[seq(1, length(endpts), by = 240)]
M.240 <- period.apply(zt, INDEX = endpts3, FUN = max)
library(qrmtools)
library(QRM)
fit.60 <- fit_GEV_MLE(M.60)
stopifnot(fit.60$convergence == 0)
xi.60 <- fit.60$par[["shape"]]
mu.60 <- fit.60$par[["loc"]]
sig.60 <- fit.60$par[["scale"]]
fit.60$SE # standard errors
fit.120 <- fit_GEV_MLE(M.120)
stopifnot(fit.120$convergence == 0)
xi.120 <- fit.120$par[["shape"]]
mu.120 <- fit.120$par[["loc"]]
sig.120 <- fit.120$par[["scale"]]
fit.120$SE # standard errors
fit.240 <- fit_GEV_MLE(M.240)
stopifnot(fit.240$convergence == 0)
xi.240 <- fit.240$par[["shape"]]
mu.240 <- fit.240$par[["loc"]]
sig.240 <- fit.240$par[["scale"]]
fit.240$SE # standard errors
# 2.2
# return level r_240,10
qGEV(1-1/10, xi.240, mu.240, sig.240)
# return period k_240,0.05
1/(1-pGEV(0.05, xi.240, mu.240, sig.240))
# 2.3
alpha <- 0.95
# VaR for Z_t
VaR.60 <- qGEV(alpha^60, xi.60, mu.60, sig.60)
VaR.120 <- qGEV(alpha^120, xi.120, mu.120, sig.120)
VaR.240 <- qGEV(alpha^240, xi.240, mu.240, sig.240)
# VaR for x_t+1
set.seed(271)
v.60 <- -0.0008441 + 0.02099*VaR.60
v.120 <- -0.0008441 + 0.02099*VaR.120
v.240 <- -0.0008441 + 0.02099*VaR.240



# Q3
# find the threshold
mean_excess_plot(zt[zt>0])
u <- 0.035
abline(v = u)
GPD_shape_plot(zt)
abline(v = u)
abline(h = 0.15)
# fit a GPD to the excesses
exceed <- zt[zt>u]
excess <- exceed - u
fit <- fit_GPD_MLE(excess)
shape.u <- fit$par[["shape"]]
scale.u <- fit$par[["scale"]]
# VaR & ES for Zt
VaR.u <- VaR_GPDtail(alpha, threshold = u, p.exceed = mean(zt>u), shape = shape.u, scale = scale.u)
ES.u <- ES_GPDtail(alpha, threshold = u, p.exceed = mean(zt>u), shape = shape.u, scale = scale.u)
# VaR & ES for x_t+1
v1 <- -0.0008441 + 0.02099*VaR.u
es1 <- -0.0008441 + 0.02099*ES.u




# Q4
par(mfrow = c(1,2))
hillPlot(zt, option = "alpha")
abline(v = 160)
abline(h = 3.1)
hillPlot(zt, option = "alpha", start = 10, end = 170)
abline(h = 3.1)
k <- 160
z.vec <- as.vector(zt[zt > 0])
alpha_hill <- 1/(mean(log(sort(z.vec, decreasing = TRUE)[1:k]))-log(sort(z.vec, decreasing = TRUE)[k]))
# VaR & ES for Zt
VaR.h <- ((2099*(1-alpha)/k)^(-1/alpha_hill))*sort(z.vec, decreasing = TRUE)[k]
ES.h <- (1+1/(alpha_hill - 1))*VaR.h
# VaR & ES for x_t+1
v2 <- -0.0008441 + 0.02099*VaR.h
es2 <- -0.0008441 + 0.02099*ES.h

