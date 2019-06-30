###################用协方差法求beta################
#步骤
#获取谷歌，标普500，一个月美元libor利率
#计算各自的收益率
library(quantmod)
getSymbols("GooG", from = '2009-06-01', to = '2013-06-01')
getSymbols("^GSPC", from = '2009-06-01', to = '2013-06-01')
getSymbols("USD1MTD156N", from = '2009-06-01', to = '2013-06-01', src = 'FRED')
LIBOR = na.omit(USD1MTD156N)
G = data.frame(GOOG$GOOG.Close)
G['date'] = index(GOOG)
sp500 = data.frame(GSPC$GSPC.Adjusted)
sp500['date'] = index(GSPC)
LIBOR = data.frame(LIBOR)
LIBOR['date'] = index(na.omit(USD1MTD156N))
#获取共同日期
sapply(list(G$date, sp500$date, LIBOR$date), length)
NG = merge(G, sp500, by = 'date')
new3 = merge(NG, LIBOR, by = 'date')
#求对数收益率
logreturn = function(x){log(tail(x, -1)/head(x, -1))}####书里是Tail/head
#一个月美元libor利率的年化收益率
#head去最后一行
rft = log(1 + head(new3$USD1MTD156N, -1)/36000*diff(as.numeric(new3$date)))
#简单beta估计
cov((logreturn(new3$GOOG.Close) - rft), (logreturn(new3$GSPC.Adjusted) - rft))/var(logreturn(new3$GSPC.Adjusted) - rft)
##############基于线性回归估计贝塔##############'
Griskpremium = logreturn(new3$GOOG.Close) - rft
Mriskpremium = logreturn(new3$GSPC.Adjusted) - rft
fit<-lm(Griskpremium ~ Mriskpremium)
fit
#0.8961131和简单beta估计得到的结果一样和书中结果0.8997有略微差异
plot(Mriskpremium, Griskpremium)
abline(fit, col = 'red')
#CAPM假设alpha为0, 模型放入-1，设置alpha为0
fit2<-lm(Griskpremium ~ -1 + Mriskpremium)#omittingintercept
summary(fit2)
summary(fit)
par(mfrow = c(2, 2))
plot(fit)
##################假设检验###################
#股票风险溢价和beta是否有显著的关系########
#使用2003年到2007年的样本换成月度收益率
library(tseries)
library(quantmod)
#已剔除时间不完整的股票
symbols<-c("A", "AA", "AAPL", "ABC", "ABT", "ACN", "ADBE", "ADI", "ADM", "ADP", 
          "ADSK", "AEE", "AEP", "AES", "AFL", "AGN", "AIG", "AIV", "AKAM", 
          "ALL", "ALXN", "AMAT", "AMD", "AMGN", "AMT", "AMZN", "AN", "ANF", 
          "AON", "APA", "APC", "APD", "APH", "ATI", "AVB", "AVP", "AVY", "AXP", 
          "AZO", "BA", "BAC", "BAX", "BBBY", "BBT", "BBY", "BDX", "BEN", "BIIB", 
          "BK", "BLK", "BMS", "BMY", "BSX", "BXP", "CAG", "CAH", "BP", "C", 
          "CAT", "CB", "CCI", "CCL", "CELG", "CERN", "CHRW", "CI", 
          "CINF", "CL", "CLX", "CMA", "CMCSA", "CME", "AMD", "AXP", "BAC", "JNJ", 
          "CSCO", "CVX", "MRK", "MSFT", "WMT", "INTC", "GS", "TRV", "IBM", "UL", 
          "VZ", "WBA", "DIS", "T", "WFC", "TM", "DIS", "UN", "E", "UNP", "CAT", "XOM", 
          "JNJ", "KO", "JPM")
res <- lapply(symbols, function(symbol)get.hist.quote(symbol, quote = "AdjClose", quiet = TRUE, start = as.Date('2003-01-01'), end = as.Date('2007-01-01')))
date = c(as.Date('2003-01-01') : as.Date('2007-01-01'))
#灰色的一列是rownames
res = as.data.frame(res)
colnames(res) = symbols
res['date'] = list(rownames(res))
library(quantmod)
getSymbols("^GSPC", from = '2003-01-01', to = '2007-01-01')
new500 = data.frame(GSPC)
new500['date'] = list(rownames(new500))
getSymbols("USD1MTD156N", from = '2003-01-01', to = '2007-01-01', src = 'FRED')
newlibor = data.frame(USD1MTD156N)
newlibor['date'] = list(rownames(newlibor))
cdate = intersect(new500$date, newlibor$date)
#提取每个月月初的时间
d = data.frame('date' = as.Date(cdate, "%Y-%m-%d"))
d$day = format(d$date, format = '%d')
d$my = format(d$date, format = '%Y-%m')
#按月份分组对date取每月最小一天
(fds = with(d, tapply(day, my, min)))
(fds = as.Date(paste(row.names(fds), fds, sep = '-')))
#取月初的数据
New500 = as.data.frame(new500[as.Date(new500$date)%in%as.Date(fds), 'GSPC.Adjusted'])
Newlibor = as.data.frame(newlibor[as.Date(newlibor$date)%in%as.Date(fds), 'USD1MTD156N'])
ress = as.data.frame(res[as.Date(res$date)%in%as.Date(fds), symbols])

#1个月美元libor有缺失值，填充
rft = as.data.frame(log(1 + head(Newlibor[, 1], -1)/36000*as.numeric(diff(fds))))
rft[17, 1] = (rft[15, 1] + rft[16, 1] + rft[18, 1] + rft[19, 1])/4
rft[25, 1] = (rft[23, 1] + rft[24, 1] + rft[26, 1] + rft[27, 1])/4
rft[29, 1] = (rft[27, 1] + rft[28, 1] + rft[30, 1] + rft[31, 1])/4
rft[41, 1] = (rft[39, 1] + rft[40, 1] + rft[42, 1] + rft[43, 1])/4
################
logreturn = function(x){log(tail(x, -1)/head(x, -1))}
new500_return = as.data.frame(apply(New500, 2, logreturn))
#求70个股票超额收益率
resss = as.data.frame(apply(ress, 2, logreturn))
a = as.data.frame(rep(rft, 101))
riskpremium = resss-a
#sp500超额收益
premium500 = new500_return-rft
#计算beta
beta = data.frame()
r = t(sapply(symbols, function(symbol)
    c(beta = lm(riskpremium[, symbol] ~ premium500[, 1])$coefficients[[2]], 
    mean = mean(resss[, symbol]))))
r = as.data.frame(r)
par(mfrow = c(1, 1))
plot(r$beta, r$mean)
abline(lm(r$mean ~ r$beta), col = 'red')
summary(lm(r$mean ~ r$beta))
#检验个体方差的解释能力
#加入非系统风险作为第二个解释变量
k = t(sapply(symbols, function(symbol){
    stock<-riskpremium[, symbol]
    beta = lm(stock ~ premium500[, 1])$coefficients[[2]]
    c(beta = beta, 
    mean = mean(stock), 
    risk = var(stock)-beta^2*var(New500)
)}))
k<-as.data.frame(k)
summary(lm(k$mean ~ k$beta + k$risk))
#股票特质风险对股票超额收益率影响虽然显著，但其值量级只有10e-7
