library(fOptions)
GBSOption(TypeFlag = "c", S = 900, X =950, Time = 1/4, r = 0.02, sigma = 0.22, b = 0.02)

GBSOption(TypeFlag = "p", S = 900, X =950, Time = 1/4, r = 0.02, sigma = 0.22, b = 0.02)@price

CRRBinomialTreeOption(TypeFlag = "ce", S = 900, X = 950, Time = 1/4, 
                      r = 0.02, b = 0.02, sigma = 0.22, n = 3)@price
CRRBinomialTreeOption(TypeFlag = "pe", S = 900, X = 950, Time = 1/4, 
                      r = 0.02, b = 0.02, sigma = 0.22, n = 3)@price

CRRTree <- BinomialTreeOption(TypeFlag = "ce", S = 900, X = 950, Time = 1/4, 
                              r = 0.02, b = 0.02, sigma = 0.22, n = 3)
BinomialTreePlot(CRRTree, dy = 1, xlab = "Time steps", ylab = "Number of up steps", 
                 xlim = c(0,4))
title(main = "Call Option Tree")

prices <- sapply(1:200, function(n) {
  CRRBinomialTreeOption(TypeFlag = "ce", S = 900, X = 950, Time = 1/4, 
                        r = 0.02, b = 0.02, sigma = 0.22, n = n)@price
})
str(prices)
price <- GBSOption(TypeFlag = "c", S = 900, X = 950, Time = 1/4, 
                   r = 0.02, sigma = 0.22, b = 0.02)@price
plot(1:200, prices, type='l', xlab = 'Number of steps', ylab = 'Prices')
abline(h = price, col ='red')
legend("bottomright", legend = c('CRR-price', 'BS-price'), col = c('black', 'red'), pch = 19)

sapply(c('delta', 'gamma', 'vega', 'theta', 'rho'), function(greek)
  GBSGreeks(Selection = greek, TypeFlag = "c", S = 900, X = 950,
            Time = 1/4, r = 0.02, b = 0.02, sigma = 0.22))

sapply(c('delta', 'gamma', 'vega', 'theta', 'rho'), function(greek)
  GBSGreeks(Selection = greek, TypeFlag = "c", S = 900, X = 950,
            Time = 1/4, r = 0.02, b = 0.02, sigma = 0.22))

deltas <- sapply(c(1/4, 1/20, 1/50), function(t)
  sapply(500:1500, function(S)
    GBSGreeks(Selection = 'delta', TypeFlag = "c",
              S = S, X = 950, Time = t, r = 0.02, b = 0.02, sigma = 0.22)))

plot(500:1500, deltas[, 1], ylab = 'Delta of call option', xlab = "Price of the underlying (S)", type = 'l')
lines(500:1500, deltas[, 2], col='blue')
lines(500:1500, deltas[, 3], col='red')
legend("bottomright", legend = c('t=1/4', 't=1/20', 't=1/50'), col = c('black', 'blue', 'red'), pch = 19)

straddles <- sapply(c('c', 'p'), function(type)
  sapply(500:1500, function(S)
    GBSGreeks(Selection = 'delta', TypeFlag = type, S = S, 
              X = 950, Time = 1/4, r = 0.02, 
              b = 0.02, sigma = 0.22)))

plot(500:1500, rowSums(straddles), type='l', xlab='Price of the underlying (S)', 
     ylab = 'Delta of straddle')

goog <- read.csv('goog_calls.csv')
volatilites <- sapply(seq_along(goog$Strike), function(i)
  GBSVolatility(price = goog$Ask.Price[i], TypeFlag = "c",
                S = 866.2, X = goog$Strike[i], Time = 88/360, 
                r = 0.02, b = 0.02))

str(volatilites)
plot(x = goog$Strike, volatilites, type = 'p', ylab = 'Implied volatiltiy', 
     xlab = 'Strike price (X)')

GBSOption(TypeFlag = "c", S = 900, X =950, Time = 1/4, 
          r = 0.02, sigma = 0.22, b = 0.02)

D2_Wiener <- function(){
  windows(200, 75)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
  for(i in 1:3){
    W1 <- cumsum(rnorm(100000))
    W2 <- cumsum(rnorm(100000))
    plot(W1,W2, type= "l", ylab = "", xlab = "")
  }
  mtext("2-dimensional Wiener-processes", outer = TRUE, 
        cex = 1.5, line = -1)
  D2_Wiener()
}

Correlated_Wiener <- function(cor){
  windows(200, 75)
  par(mfrow = c(1, 3), oma = c(0, 0, 2, 0))
  for(i in 1:3){
    W1 <- cumsum(rnorm(100000))
    W2 <- cumsum(rnorm(100000))
    W3 <- cor*W1+sqrt(1-cor^2)*W2
    plot(W1,W3, type= "l", ylab = "", xlab = "")
  }
  mtext("Correlated Wiener-processes", outer = TRUE, cex = 1.5, line = -1)
}
Correlated_Wiener(0.6)

Margrabe <- function(S1, S2, sigma1, sigma2, Time, rho, delta1 = 0, delta2 = 0){
  sigma <- sqrt(sigma1^2 + sigma2^2 - 2*sigma1*sigma2*rho)
  d1 <- ( log(S1/S2) + ( delta2-delta1 + sigma^2/2 ) * Time ) / (sigma*sqrt(Time))
  d2 <- ( log(S1/S2) + ( delta2-delta1 - sigma^2/2 ) * Time ) / (sigma*sqrt(Time))
  M <- S1*exp(-delta1*Time)*pnorm(d1) - S2*exp(-delta2*Time)*pnorm(d2)
  return(M)
}
Margrabe(100, 120, 0.2, 0.3, 2, 0.15)

x <- seq(-1, 1, length = 1000)
y <- rep(0, 1000)
for(i in 1:1000)
  y[i] <- Margrabe(100, 120, 0.2, 0.3, 2, x[i])
plot(x, y, xlab = "correlation", ylab = "prices",
     main = "Price of exchange option", type = "l", lwd = 3)        

vasicek = function(alpha, beta, sigma, n = 1000, r0 = 0.05){
  v <- rep(0, n)
  v[1] = r0
  for (i in 2:n){
    v[i] <- v[i-1]+alpha*(beta - v[i-1]) + sigma*rnorm(1)
  }
  return(v)
}


r = matrix(0, 1000, 3)
set.seed(2014)
r[,1] <- vasicek(0.002, 0.065, 0.0002)
set.seed(2014)
r[,2] <- vasicek(0.02, 0.065, 0.0002)
set.seed(2014)
r[,3] <- vasicek(0.2, 0.065, 0.0002)

matplot(r, type = "l", ylab = "", xaxt = "no",  main = "Vasicek trajectories with alpha = 0.2%, 2% and 20%")
lines(c(-1,1001), c(0.065, 0.065), col = "grey", lwd = 2, lty = 3)

vasicek_pdf = function(x, alpha, beta, sigma, delta_T, r0 = 0.05){
  e <- r0*exp(-alpha*delta_T)+beta*(1-exp(-alpha*delta_T))
  s <- sigma^2/(2*alpha)*(1-exp(-2*alpha*delta_T))
  y <- dnorm(x, mean = e, sd = s)
  return(y)
}

x <- seq(-0.1, 0.2, length = 1000)
y1 <- vasicek_pdf(x, .2, 0.1, 0.15,10)
y2 <- vasicek_pdf(x, .2, 0.1, 0.15, 5)
y3 <- vasicek_pdf(x, .2, 0.1, 0.15, 3)
y4 <- vasicek_pdf(x, .2, 0.1, 0.15, 2)

par(xpd = T ,mar = c(2,2,2,2), mfrow = c(2,2))  
matplot(x, cbind(y1,y2,y3,y4), type = "l",ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topleft", c("T-t = 2", "T-t = 3", "T-t = 5", "T-t = 10"), fill = 1:5, cex = 0.7)

y1 <- vasicek_pdf(x, .2, 0.1, 0.15, 5)
y2 <- vasicek_pdf(x, .2, 0.12, 0.15, 5)
y3 <- vasicek_pdf(x, .2, 0.14, 0.15, 5)
y4 <- vasicek_pdf(x, .2, 0.16, 0.15, 5)

matplot(x, cbind(y1,y2,y3,y4), type = "l", ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topleft", c("beta = 0.1", "beta = 0.12", "beta = 0.14", "beta = 0.16"), fill = 1:5, cex = 0.7)

y1 <- vasicek_pdf(x, .1, 0.1, 0.15, 5)
y2 <- vasicek_pdf(x, .2, 0.1, 0.15, 5)
y3 <- vasicek_pdf(x, .3, 0.1, 0.15, 5)
y4 <- vasicek_pdf(x, .4, 0.1, 0.15, 5)

matplot(x, cbind(y1,y2,y3,y4), type = "l", ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topleft", c("alpha = 0.1", "alpha = 0.2", "alpha = 0.3", "alpha = 0.4"), fill = 1:5, cex = 0.7)

y1 <- vasicek_pdf(x, .1, 0.1, 0.10, 5)
y2 <- vasicek_pdf(x, .1, 0.1, 0.12, 5)
y3 <- vasicek_pdf(x, .1, 0.1, 0.14, 5)
y4 <- vasicek_pdf(x, .1, 0.1, 0.15, 5)

matplot(x, cbind(y1,y2,y3,y4), type = "l", ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topleft", c("sigma = 0.1", "sigma = 0.12", "sigma = 0.14", "sigma = 0.15"), fill = 1:5, cex = 0.7)

CIR_pdf = function(x, alpha, beta, sigma, delta_T, r0 = 0.1){
  q = (2*alpha*beta)/(sigma^2) - 1
  c = (2*alpha)/(sigma^2*(1-exp(-alpha*delta_T)))
  u = c*r0*exp(-alpha*delta_T)
  y = dchisq(2*c*x, df = 2*q+2, ncp = 2*u)
}


x <- seq(0, 0.15, length = 1000)
y1 <- CIR_pdf(x, .3, 0.05, 0.1, 1)
y2 <- CIR_pdf(x, .3, 0.05, 0.1, 2)
y3 <- CIR_pdf(x, .3, 0.05, 0.1, 5)
y4 <- CIR_pdf(x, .3, 0.05, 0.1, 50)

par(mar = c(2,2,2,2), mfrow = c(2,2))  
matplot(x, cbind(y1,y2,y3,y4), type = "l",ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topright", c("T-t = 1", "T-t = 2", "T-t = 5", "T-t = 50"), fill = 1:5, cex = 0.7)

y1 <- CIR_pdf(x, .2, 0.05, 0.1, 1)
y2 <- CIR_pdf(x, .4, 0.05, 0.1, 1)
y3 <- CIR_pdf(x, .6, 0.05, 0.1, 1)
y4 <- CIR_pdf(x,  1, 0.05, 0.1, 1)

matplot(x, cbind(y1,y2,y3,y4), type = "l",ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topright", c("alpha = 0.2", "alpha = 0.4", "alpha = 0.6", "alpha = 1"), fill = 1:5, cex = 0.7)

y1 <- CIR_pdf(x, .3, 0.10, 0.1, 1)
y2 <- CIR_pdf(x, .3, 0.12, 0.1, 1)
y3 <- CIR_pdf(x, .3, 0.14, 0.1, 1)
y4 <- CIR_pdf(x, .3, 0.16, 0.1, 1)

matplot(x, cbind(y1,y2,y3,y4), type = "l",ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topright", c("beta = 0.1", "beta = 0.12", "beta = 0.14", "beta = 0.16"), fill = 1:5, cex = 0.7)

x <- seq(0, 0.25, length = 1000)
y1 <- CIR_pdf(x, .3, 0.05, 0.03, 1)
y2 <- CIR_pdf(x, .3, 0.05, 0.05, 1)
y3 <- CIR_pdf(x, .3, 0.05, 0.10, 1)
y4 <- CIR_pdf(x, .3, 0.05, 0.15, 1)

matplot(x, cbind(y1,y2,y3,y4), type = "l",ylab ="",xlab = "", col = 1:5, lty = 1)
legend("topright", c("sigma = 3%", "sigma = 5%", "sigma = 10%", "sigma = 15%"), fill = 1:5, cex = 0.7)

library(SMFI5)
bond.vasicek(0.5,2.55,0.365,0.3,0,3.55,1080,c(1/12,3/12,6/12,1),365)