library(TSA)

data <- read.csv("papers_published_per_month.csv", header = T, row.names = 1)

# plot the time series
complete_ts <- ts(data$count, frequency = 12, start = min(data$year))
plot(complete_ts, ylab = "Number of Papers Published", type = "l", xlab = "Year")

# Zoom in on 2005-2010
partial_ts <- ts(data[which(data$year>=2005 & data$year<=2010),3], frequency = 12, start =2005)
plot(partial_ts, ylab = "Number of Papers Published", type = "l", xlab = "Year")
points (y = partial_ts, x = time(partial_ts), pch= as.vector(season(partial_ts)), cex = 1.5)

# plot stl
plot(stl(complete_ts, "periodic"), xlab = 'Year')

# estimate order of Box Cox transformation
BoxCox.ar(complete_ts,method = 'ols',lambda = seq(0.36,.65,0.01))

# plot first difference
first_diff_ts <- diff(sqrt(complete_ts))
plot (first_diff_ts, ylab = "First Difference of Square Root Counts", xlab = "Year")
abline(h = 0, col = "red")
par(mfrow = c(1,2))
acf (first_diff_ts, lag.max = 50, main = "")
pacf (first_diff_ts,lag.max = 50, main = "")

# plot first and seasonal differences
first_and_seasonal_diff_ts <- diff(diff(sqrt(complete_ts)),lag = 12)
par(mfrow = c(1,3))
years = rep (seq(1,12,1), 48)
boxplot (first_diff_ts ~ years[-1] , xlab = "Month", ylab = "First Difference of Square Root", main = "Before Seasonal Difference")
boxplot(first_and_seasonal_diff_ts ~ years[-c(1:13)], xlab = "Month", ylab = "First and Seasonal Difference of Square Root", main = "After Seasonal Difference")
plot (first_and_seasonal_diff_ts, ylab = "First and Seasonal Difference of Square Root Counts", xlab = "Year", main = "First and Seasonal Difference")
abline( h = 0, col = "red")
par(mfrow = c(1,2))
acf(first_and_seasonal_diff_ts, lag.max = 40, main = "")
pacf(first_and_seasonal_diff_ts,lag.max = 60, main = "")

# determine optimal model parameters
coefs <- expand.grid(p = seq(0,3),P = seq(0,3),q = seq(0,3),Q = seq(0,3))
aic_values <- data.frame()
for (i in 1:dim(coefs)[1]){
  mod = try(arima(sqrt(complete_ts), order = c(coefs$p[i],1,coefs$q[i]),seasonal=list(order=c(coefs$P[i],1,coefs$Q[i]),period=12)))
  if (is(mod,"try-error")){
    print(paste(i, "- Error fitting model for parameters:", coefs$p[i], coefs$P[i], coefs$q[i], coefs$Q[i]))
  }
  else {
    aic <- mod$aic
    aic_values <- rbind(aic_values, cbind(coefs[i,], aic))
    print(paste(i, "- AIC calculated for parameters:", coefs$p[i], coefs$P[i], coefs$q[i], coefs$Q[i]))
  }
}
head(aic_values[order(aic_values$aic),])

# fit model
arima_model = arima(sqrt(complete_ts), order = c(2,1,1), seasonal=list(order=c(2,1,3),period=12))
detectIO(arima_model)
final_model = arima(sqrt(complete_ts), order=c(2,1,1), seasonal=list(order=c(2,1,3),period=12), io = c(37, 565))

# analysis of residuals
TSA::spec(final_model$residuals, main = "Periodogram of Residuals", plot=T)
par(mfrow=c(1,2))
plot(final_model$residuals ~ fitted(final_model), main='Residuals Vs. Fitted', ylab = "Residual", xlab = "Fitted Value")
abline( h = 3, lty= 'dotted', col = 'blue')
abline (h = -3, lty= 'dotted', col = 'blue')
abline (h = 0, col = 'red')
acf(final_model$residuals, main='Sample ACF of Residuals', lag.max=20)
qqnorm(final_model$residuals, main='Q-Q Plot of Estimated Residuals')
qqline(final_model$residuals, col = "red")

# Ljung-Box test
pval <- vector(length = 20)
for (i in 1:length(pval)){ 
  b1 <- Box.test(final_model$residuals, lag=i, type="Ljung") 
  pval[i]<- b1$p.value
}
pval<0.05

# spectral analysis
sp = TSA::spec(first_and_seasonal_diff_ts, spans = 8, lty = 'dashed', ci.plot=T, col = "blue", main = "")
f = seq(0.001,0.5,by = 0.009)
lines (f, ARMAspec(model = list(ar = final_model$coef[grep("^ar", names(final_model$coef))], ma = final_model$coef[grep("^ma", names(final_model$coef))], seasonal = list(sar=  final_model$coef[grep("^sar", names(final_model$coef))], sma = final_model$coef[grep("^sma", names(final_model$coef))], period = 12)), freq = f, plot = F)$spec, lty="solid", col = "red")
