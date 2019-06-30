library(ggplot2)
library(moments)
library(rugarch)
library(zoo)
library(qrmtools)

return <- read.csv("同仁堂收益率.csv", header = TRUE)
return$日期 <- as.Date(return$日期)
ggplot(return, aes(x = 日期, y = 收益率)) + geom_line() + theme_minimal()

# 偏度、峰度计算
skewness(return$收益率)
kurtosis(return$收益率)
ggplot(return, aes(收益率, fill = I('lightblue'))) + geom_histogram(bins=100) + theme_minimal()

# 正态性检验
jarque.test(return$收益率)
qqnorm(return$收益率)
qqline(return$收益率)

# 自相关图
acf(return$收益率,main="收益率自相关图")
acf(abs(return$收益率), main="收益率绝对值自相关图")

# Ljung-Box检验
for(i in seq(2,10,2)){
  print(Box.test(return$收益率, lag = i, type="Ljung-Box"))
}
for(i in seq(2,10,2)){
  print(Box.test(abs(return$收益率), lag = i, type="Ljung-Box"))
}

# 高斯新息模型
GARCHnorm.spec <- ugarchspec(mean.model = list(armaOrder = c(3,3), include.mean = FALSE),
                             variance.model = list(garchOrder = c(1,1)),
                             distribution.model = "norm")
model.norm <- ugarchfit(GARCHnorm.spec, return$收益率)
Box.test(residuals(model.norm, standardize=TRUE)^2, lag = 6, type="Ljung-Box")

layout(matrix(1:4,ncol = 2, nrow = 2, byrow = TRUE))
plot(model.norm, which = 9) # Q-Q图
plot(model.norm, which = 10) # ACF图
# t分布新息模型
GARCHt.spec <- ugarchspec(mean.model = list(armaOrder = c(3,3), include.mean = FALSE),
                             variance.model = list(garchOrder = c(1,1)),
                             distribution.model = "std")
model.t <- ugarchfit(GARCHt.spec, return$收益率)
Box.test(residuals(model.t, standardize=TRUE)^2, lag = 6, type="Ljung-Box")
plot(model.t, which = 9) # Q-Q图
plot(model.t, which = 10) # ACF图

# 模型比较
(AIC.norm <- 2 * length(model.norm@fit$coef) - 2 * model.norm@fit$LLH)
(AIC.t <- 2 * length(model.t@fit$coef) - 2 * model.t@fit$LLH)

(BIC.norm <- log(length(return$收益率)) * length(model.norm@fit$coef) - 2 * model.norm@fit$LLH)
(BIC.t <- log(length(return$收益率)) * length(model.t@fit$coef) - 2 * model.t@fit$LLH)

# 预测
forc.norm <- ugarchforecast(model.norm, n.head = 1)
fitted(forc.norm)[1]
sigma(forc.norm)[1]
as.numeric(quantile(forc.norm, 0.95))[1] # VaR
ES_t(0.95, scale = sigma(forc.norm), df = Inf) # ES

forc.t <- ugarchforecast(model.t, n.head = 1)
fitted(forc.t)[1]
sigma(forc.t)[1]
as.numeric(quantile(forc.t, 0.95))[1] # VaR
ES_t(0.95, scale = sigma(forc.norm), df = model.t@fit$coef["shape"]) # ES

