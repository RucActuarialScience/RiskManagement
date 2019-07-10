library(RQuantLib)
# 使用的的日期定义为今天
today <- Sys.Date() 
# 设置交易,结算日期和远期利率时间间隔
params <- list(tradeDate = today - 2, settleDate = today, dt = 0.25)
# 给定一条平坦的收益率曲线,从0 到10 步长为0.1 的序列
times <- seq(0, 10, 0.1) 
# 股息率为0
dividendYield <- DiscountCurve(params, list(flat = 10e-6), times)
# 无风险收益率为5%
riskFreeRate <- DiscountCurve(params, list(flat = 0.05), times)

#构建BlackScholes过程和二叉树定价的参数
process <- list(underlying = 20, divYield = dividendYield,
                rff = riskFreeRate, volatility = 0.2)

#构建可转换证券部分的参数，包括期权类型、证券面值、可赎回价格、利差、转换
bondparams <- list(exercise = "eu", faceAmount = 100, redemption = 100, creditSpread = 0.02,
                   conversionRatio = 4, issueDate = as.Date(today + 2),
                   maturityDate = as.Date(today + 1825))
#我们还需要设定转换比率，它决定了如果债券持有者决定把债券转换为股票，那他会得到多少普通股。同时也在这里设定债券的票面价值和信用利差。
dateparams <- list(settlementDays = 3, dayCounter = "ActualActual",
                   period = "Annual", businessDayConvention = "Unadjusted")
ConvertibleFixedCouponBond(bondparams, coupon = 0.05, process, dateparams)

# 现在，我们从1～30 提高基础股票的价格，同时检查净现值的改变。
res <- sapply(seq(1, 30, 1), function(s) {
  process$underlying = s
  ConvertibleFixedCouponBond(bondparams, coupon = 0.05, process,
                             dateparams)$NPV
})
plot(1:30, res, type = 'l', xlab = 'Price of the underlying stocks', 
     ylab = 'Net Present Value')