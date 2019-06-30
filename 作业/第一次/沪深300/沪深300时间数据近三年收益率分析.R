library(rugarch)
library(zoo)
library(ADGofTest) # for ad.test()
library(moments) # for skewness(), kurtosis()
library(qrmtools)
library(xts)
library(ggplot2)

### input data #############################################################
mydata<-read_excel("data.xlsx")
dat1=as.matrix(mydata[3410:4145,3])
dat1=as.numeric(dat1)
dat1=dat1/100
Time=as.matrix(mydata[3410:4145,1])
mm=mydata[3410:4145,-2]
mm=mm[,-3]
names(mm)<-c("Time","rate")
rate <- xts(mm$rate, as.Date(mm$Time, format='%Y/%m/%d'))
rate=rate/100
png(file='F:/risk_figure/rate.png')
plot.zoo(rate, xlab = "Time")
dev.off()
############compute ##########

piandu=skewness(dat1)
fengdu=kurtosis(dat1)

############ QQ figure #################
alpha=0.05 #alpha为显著性水平，这里的默认值为0.05 
pdf(file='F:/risk_figure/qq_dat1.pdf')
# par(mfrow=c(2,1))
qqnorm(dat1,main="Q-Q")
qqline(dat1)
hist(dat1,freq=F,main="Histograms and density estimation curves")
lines(density(dat1),col="blue")#密度估计曲线
x<-c(round(min(dat1)):round(max(dat1)))
lines(x,dnorm(x,mean(dat1),sd(dat1)),col="red")#正态分布曲线
dev.off()
sol<-shapiro.test(dat1)
if(sol$p.value>alpha){
  print(paste("success:服从正态分布,p.value=",sol$p.value,">",alpha))
}else{
  print(paste("error:不服从正态分布,p.value=",sol$p.value,"<=",alpha))
}
sol
############## Jarque-Bera检验 ############
(sh <- shapiro.test(dat1)) # Shapiro--Wilk
(ag <- agostino.test(dat1)) # D'Agostino's test
(jb <- jarque.test(dat1)) # Jarque--Bera test

###############ACF图像#################

bacf <- acf(dat1, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

png(file='F:/risk_figure/acf_dat1.png')
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_segment(mapping = aes(xend = lag, yend = 0),color='blue',size=5,alpha=I(1/2)) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, color = 'darkblue')+
  geom_hline(aes(yintercept=0))
dev.off()


#abs dat1


abs_dat1=abs(dat1)
bacf <- acf(abs_dat1, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

png(file='F:/risk_figure/acf_abs_dat1.png')
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  geom_segment(mapping = aes(xend = lag, yend = 0),color='blue',size=5,alpha=I(1/2)) +
  geom_hline(aes(yintercept = 0.05), linetype = 2, color = 'darkblue')+
  geom_hline(aes(yintercept=0))
dev.off()


###############Ljung-Box #################
LB.raw   <- apply(rate,        2, Box.test, lag = 10, type = "Ljung-Box")
LB.abs   <- apply(abs(rate),   2, Box.test, lag = 10, type = "Ljung-Box") 

### 2 Fit an ARMA(1,1)--GARCH(1,1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,1), # AR(1,1) part
                                        include.mean = TRUE), # with mean
                      distribution.model = "norm") # normal innovations
(fit.N <- ugarchfit(spec = uspec.N, data = rate))

png(file='F:/risk_figure/model_nn.png')
layout(matrix(1:2, ncol = 1, byrow = TRUE)) # specify (2,2)-matrix of plots
plot(fit.N, which = 9) # Q-Q plot of standardized residuals Z_t (shows leptokurtosis of Z_t; normal assumption not supported)
plot(fit.N, which = 10) # ACF of standardized residuals Z_t (shows AR dynamics do a reasonable job of explaining conditional mean)
layout(1) # restore layout
dev.off()

### 3 Fit an AR(1,1)--GARCH(1,1) model with Student t innovations ################

## Now consider t innovations
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                      mean.model = list(armaOrder = c(1,1), include.mean = TRUE),
                      distribution.model = "std") # Student t innovations
fit.t <- ugarchfit(spec = uspec.t, data = rate)

## The pictures are similar, but the Q-Q plot looks "better"
png(file='F:/risk_figure/model_t.png')
layout(matrix(1:2, ncol = 1, byrow = TRUE)) # specify (2,2)-matrix of plots
plot(fit.t, which = 9) # Q-Q plot of Z_t against a normal
plot(fit.t, which = 10) # ACF of Z_t
layout(1)
dev.off()

###  Using AIC BIC to compare the two models ###########

LL.N <- fit.N@fit$LLH # log-likelihood of the model based on normal innovations
LL.t <- fit.t@fit$LLH
(AIC.N <- 2*length(fit.N@fit$coef) - 2*LL.N)
(AIC.t <- 2*length(fit.t@fit$coef) - 2*LL.t)
## => t model is preferred
(BIC.N <- log(736)*length(fit.N@fit$coef) - 2*LL.N)
(BIC.t <- log(736)*length(fit.t@fit$coef) - 2*LL.t)


### Forecasting/predicting from the fitted model #############################
forc_n <- ugarchforecast(fit.N, n.ahead = 1)
forc_t <- ugarchforecast(fit.t, n.ahead = 1)

#VaR and ES fit.t

param.t <- coef(fit.t) # estimated coefficients
sig.t <- sigma(fit.t) # estimated volatility
VaR.95 <- quantile(fit.t, probs = 0.95) # estimated VaR at level 99%
Z.t <- residuals(fit.t, standardize = TRUE) # estimated standardized residuals Z_t

## Plots
png(file='F:/risk_figure/sig_t.png')
plot.zoo(sig.t,    xlab = "Time", ylab = expression(hat(sigma)[t]))
dev.off()
png(file='F:/risk_figure/vaR_t.png')
plot.zoo(VaR.95, xlab = "Time", ylab = expression(widehat(VaR)[0.95]))
dev.off()
png(file='F:/risk_figure/res_t.png')
plot.zoo(Z.t,      xlab = "Time", ylab = expression(hat(Z)[t]))
dev.off()

Z_alpha.t=quantile(Z.t,probs=0.95)
VaR.t_new.95 = Z_alpha.t*0.01009+0.0004877
ES.t_new.95 =  sum(Z.t*(Z.t>=Z_alpha.t))/(sum(Z.t>=Z_alpha.t))*0.01009+0.0004877

#VaR and ES fit.n

param.n <- coef(fit.N) # estimated coefficients
sig.N <- sigma(fit.N) # estimated volatility
VaR.95.N <- quantile(fit.N, probs = 0.95) # estimated VaR at level 99%
Z.N <- residuals(fit.N, standardize = TRUE) # estimated standardized residuals Z_t

## Plots
png(file='F:/risk_figure/sig_N.png')
plot.zoo(sig.N,    xlab = "Time", ylab = expression(hat(sigma)[t]))
dev.off()
png(file='F:/risk_figure/vaR_N.png')
plot.zoo(VaR.95.N, xlab = "Time", ylab = expression(widehat(VaR)[0.95]))
dev.off()
png(file='F:/risk_figure/res_N.png')
plot.zoo(Z.N,      xlab = "Time", ylab = expression(hat(Z)[t]))
dev.off()

Z_alpha.N=quantile(Z.N,probs=0.95)
VaR.N_new.95 = Z_alpha.N*0.01025-3.134e-05
ES.N_new.95 =  sum(Z.N*(Z.N>=Z_alpha.N))/(sum(Z.N>=Z_alpha.N))*0.01025-3.134e-05
