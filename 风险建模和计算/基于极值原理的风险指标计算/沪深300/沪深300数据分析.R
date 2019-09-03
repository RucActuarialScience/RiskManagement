
library(qrmdata)
library(qrmtools)
library(rugarch)
library(zoo)
library(openxlsx)
library(xts)
library(QRM)

data = read.xlsx("沪深300时间数据.xlsx", sheet = 1, detectDates = TRUE)
date = data$X1
series = data$close
hushen = zoo(series, date)
X <- -returns(hushen)
plot.zoo(X, xlab = "Time", ylab = expression(X_t))

## Question 1
### Fit an ARMA(1, 1)--GARCH(1, 1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
           mean.model = list(armaOrder = c(1, 1), # ARMA(1, 1) part
                    include.mean = TRUE), # with mean
           distribution.model = "norm") # normal innovations
(fit.N <- ugarchfit(spec = uspec.N, data = X[-1]))

### Fit an ARMA(1, 2)--GARCH(1, 1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N_add <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
           mean.model = list(armaOrder = c(1, 2), # ARMA(1, 1) part
                    include.mean = TRUE), # with mean
           distribution.model = "norm") # normal innovations
(fit.N_add <- ugarchfit(spec = uspec.N_add, data = X[-1]))

### Fit an ARMA(2, 1)--GARCH(1, 1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N_add <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
             mean.model = list(armaOrder = c(2, 1), # ARMA(1, 1) part
                      include.mean = TRUE), # with mean
             distribution.model = "norm") # normal innovations
(fit.N_add <- ugarchfit(spec = uspec.N_add, data = X[-1]))

### Fit an ARMA(1, 1)--GARCH(2, 1) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N_add <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 1)), 
             mean.model = list(armaOrder = c(1, 1), # ARMA(1, 1) part
                      include.mean = TRUE), # with mean
             distribution.model = "norm") # normal innovations
(fit.N_add <- ugarchfit(spec = uspec.N_add, data = X[-1]))

### Fit an ARMA(1, 1)--GARCH(1, 2) model with normal innovations ###################

## Model specification (without fixed.pars, so without keeping any parameter
## fixed during the optimization)
uspec.N_add <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), 
             mean.model = list(armaOrder = c(1, 1), # ARMA(1, 1) part
                      include.mean = TRUE), # with mean
             distribution.model = "norm") # normal innovations
(fit.N_add <- ugarchfit(spec = uspec.N_add, data = X[-1]))

## Series with conditional quantiles
plot(fit.N, which = 2)
layout(matrix(1:4, ncol = 2, byrow = TRUE)) # specify (2, 2)-matrix of plots
plot(fit.N, which = 6) # ACF of absolute data |X_t| (shows serial correlation)
plot(fit.N, which = 9) # Q-Q plot of standardized residuals Z_t (shows leptokurtosis of Z_t; normal assumption not supported)
plot(fit.N, which = 10) # ACF of standardized residuals Z_t (shows ARMA dynamics do a reasonable job of explaining conditional mean)
plot(fit.N, which = 11) # ACF of squared standardized residuals Z_t^2 (shows GARCH dynamics do a reasonable job of explaining conditional sd)
layout(1) # restore layout

## Now consider t innovations
uspec.t <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
           mean.model = list(armaOrder = c(1, 1), include.mean = TRUE), 
           distribution.model = "std") # Student t innovations
fit.t <- ugarchfit(spec = uspec.t, data = X[-1])

## The pictures are similar, but the Q-Q plot looks "better"
layout(matrix(1:4, ncol = 2, byrow = TRUE)) # specify (2, 2)-matrix of plots
plot(fit.t, which = 6) # ACF of |X_t|
plot(fit.t, which = 9) # Q-Q plot of Z_t against a normal
plot(fit.t, which = 10) # ACF of Z_t
plot(fit.t, which = 11) # ACF of Z_t^2
layout(1)

### Using AIC and BIC to compare the two models ###########

## Akaike's information criterion (AIC) can be used for model comparison.
## Favour the model with smallest AIC.
## Note: AIC = 2*<number of estimated parameters> - 2 * log-likelihood
LL.N <- fit.N@fit$LLH # log-likelihood of the model based on normal innovations
LL.t <- fit.t@fit$LLH # log-likelihood of the model based on t innovations
(AIC.N <- 2*length(fit.N@fit$coef) - 2*LL.N)
(AIC.t <- 2*length(fit.t@fit$coef) - 2*LL.t)
## = > t model is preferred

## BIC can be used for model comparison.
## Favour the model with smallest BIC.
## Note: BIC = log(n)*<number of estimated parameters> - 2 * log-likelihood
(BIC.N <- log(length(X))*length(fit.N@fit$coef) - 2*LL.N)
(BIC.t <- log(length(X))*length(fit.t@fit$coef) - 2*LL.t)
## = > t model is preferred

## Question 2
Z <- residuals(fit.t, standardize = FALSE) # estimated residuals Z_t
###Block Maxima Method (BMM) ################################################
#n = 60 maxima
endpts60 = c(seq(from = 0, to = length(Z), by = 60), length(Z))
Z.60 <- period.apply(Z, INDEX = endpts60, FUN = max) 
fit.60 <- fit_GEV_MLE(Z.60) # maximum likelihood estimator
stopifnot(fit.60$convergence == 0) # = > converged
(xi.60 <- fit.60$par[["shape"]]) # ~ = 0.1741
(mu.60 <- fit.60$par[["loc"]])
(sig.60 <- fit.60$par[["scale"]])
fit.60$SE # standard errors

#n = 120 maxima
endpts120 = c(seq(from = 0, to = length(Z), by = 120), length(Z))
Z.120 <- period.apply(Z, INDEX = endpts120, FUN = max) 
fit.120 <- fit_GEV_MLE(Z.120) # maximum likelihood estimator
stopifnot(fit.120$convergence == 0) # = > converged
(xi.120 <- fit.120$par[["shape"]]) # ~ = 0.0353
(mu.120 <- fit.120$par[["loc"]])
(sig.120 <- fit.120$par[["scale"]])
fit.120$SE # standard errors

#n = 240 maxima
endpts240 = c(seq(from = 0, to = length(Z), by = 240), length(Z))
Z.240 <- period.apply(Z, INDEX = endpts240, FUN = max) 
fit.240 <- fit_GEV_MLE(Z.240) # maximum likelihood estimator
stopifnot(fit.240$convergence == 0) # = > converged
(xi.240 <- fit.240$par[["shape"]]) # ~ = -0.662
(mu.240 <- fit.240$par[["loc"]])
(sig.240 <- fit.240$par[["scale"]])
fit.240$SE # standard errors

#return level
qrmtools::qGEV(1-1/10, shape = xi.240, loc = mu.240, scale = sig.240) # r_{n = 240, k = 10} ~ = 0.0810
#return period
k <- 1/(1-qrmtools::pGEV(0.05, shape = xi.240, loc = mu.240, scale = sig.240)) # corresponding estimated return period k
#VaR and ES for Z and X by empirical quantile
alpha = 0.95
(forc <- ugarchforecast(fit.t, n.ahead = 1)) #mu_(t+1) = -0.0002961, sigma_(t+1) = 0.01091
#VaR.Z = quantile(Z, probs = 0.95)
#VaR_new_emp.95 = -0.0002961+VaR.Z*0.01091
#ES_new_emp.95 =  sum(Z*(Z>VaR.Z))/(sum(Z>VaR.Z))*0.01091-0.0002961

#VaR estimation based on GEV
p = 0.05
VaR_GEV = mu.120-sig.120/xi.120*(1-(-120*log(1-p))^(-xi.120))
VaR_new_GEV = -0.0002961+VaR_GEV*0.01091

## Question 3
## Plot 
plot.zoo(Z, xlab = "time", ylab = expression(Z))
Z_excess = Z[Z>0]
par(mfrow = c(1, 2))
## Sample mean excess plot 
mean_excess_plot(Z_excess) # sample mean excess plot
u <- 0.03 # threshold 
## Effect of changing the threshold on xi
GPD_shape_plot(Z)
abline(v = 0.03)

## Fit GPD models to excesses via MLE
## u = 0.03
exceed.u <- Z[Z > u] # exceedances
excess.u <- exceed.u - u # excesses
(fit.u <- fit_GPD_MLE(excess.u)) # MLE
shape.u <- fit.u$par[["shape"]]
scale.u <- fit.u$par[["scale"]]

## Plot empiricial excess df vs fitted GPD G_{xi, beta}
## u = 0.03
res <- edf_plot(excess.u)
z.res <- tail(res$t, n = -1)
lines(z.res, col = 'red', qrmtools::pGPD(z.res, shape = shape.u, scale = scale.u)) # fitted GPD

## VaR_alpha, ES_alpha for alpha and threshold
alpha <- 0.95
(VaR.u <- VaR_GPDtail(alpha, threshold = u, p.exceed = mean(Z > u), 
            shape = shape.u, scale = scale.u))
(ES.u <- ES_GPDtail(alpha, threshold = u, p.exceed = mean(Z > u), 
           shape = shape.u, scale = scale.u))
VaR_new_GPD.95 = -0.0002961+0.01091*VaR.u
ES_new_GPD.95 = -0.0002961+0.01091*ES.u

##Hill estimator
par(mfrow = c(1, 2))
hillPlot(Z, option = "alpha")
abline(v = 150)
abline(h = 2.5)
hillPlot(Z, option = "alpha", start = 15, end = 250)
abline(h = 2.5)
k = 150
Z.vec = as.vector(Z)
alpha_hill = 1/(mean(log(sort(Z.vec, decreasing = TRUE)[1:k]))-log(sort(Z.vec, decreasing = TRUE)[k]))

#VaR and ES for hill estimator
VaR_hill.95 = (length(Z)/k*(1-alpha))^(-1/alpha_hill)*sort(Z.vec, decreasing = TRUE)[k]
ES_hill.95 = alpha_hill/(alpha_hill-1)*VaR_hill.95
VaR_new_hill.95 = -0.0002961+0.01091*VaR_hill.95
ES_new_hill.95 = -0.0002961+0.01091*ES_hill.95
