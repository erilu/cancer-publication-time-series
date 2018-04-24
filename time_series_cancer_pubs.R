# Erick Lu
# publication_ts_model.R

setwd ( "~/Documents/github/papertracker")
library(TSA)

#list of dates corresponding to pubmed abstracts
dates = read.table ("cancer_dates.txt", header = F, sep = ",")

head(dates)
years = matrix(c(1,1), nrow = 1, ncol = 2)
bymonth = matrix(c(1,1), nrow = 1, ncol = 2)

# convert dates into frequency counts, # of papers published that month of that year.
for (i in 1:58) {
  year = 1955 + (i-1)
  year12 = dates[ which (dates$V1 == year),]
  for (i in 1:12) {
    i = 12-(i-1)
    if (i == 1) {
      months = year12 [ which (year12$V2 == 'Jan'),]
    }
    if (i == 2) {
      months = year12 [ which (year12$V2 == 'Feb'),]
    }
    if (i == 3) {
      months = year12 [ which (year12$V2 == 'Mar'),]
    }
    if (i == 4) {
      months = year12 [ which (year12$V2 == 'Apr'),]
    }
    if (i == 5) {
      months = year12 [ which (year12$V2 == 'May'),]
    }
    if (i == 6) {
      months = year12 [ which (year12$V2 == 'Jun'),]
    }
    if (i == 7) {
      months = year12 [ which (year12$V2 == 'Jul'),]
    }
    if (i == 8) {
      months = year12 [ which (year12$V2 == 'Aug'),]
    }
    if (i == 9) {
      months = year12 [ which (year12$V2 == 'Sep'),]
    }
    if (i == 10) {
      months = year12 [ which (year12$V2 == 'Oct'),]
    }
    if (i == 11) {
      months = year12 [ which (year12$V2 == 'Nov'),]
    }
    if (i == 12) {
      months = year12 [ which (year12$V2 == 'Dec'),]
    }
    frequency = nrow(months)
    bymonth = rbind ( bymonth, c(year, frequency))
  }
  years = rbind(years, year12)
}

bymonth = bymonth[-1,]

files = file("cancer_per_month.txt")
write.table(bymonth, sep = ",",files)
close(files)


#Read in the data from a txt file.
data = read.table ("cancerbymonth.txt", sep = ",", header = T)

data = data[-c(685:696),]
data

data.new = data[-c(1:108),]

data.new2 = data[c(601:672),]

data.ts2 = ts (data.new2$V2, frequency = 12, start = 2005)
data.ts2


data.ts.new = ts(data.new$V2, frequency = 12, start = 1964)

#Plot the raw data with monthly symbol
par(mfrow = c(2,1))
plot(data.ts2, ylab = "Total number of Papers Published", type = "l", xlab = "Year")
points (y = data.ts2, x = time(data.ts2), pch= as.vector(season(data.ts2)), cex = 1.5)


#Plot the raw data.
par(mfrow = c(1,1))
plot(data.ts.new, main = "Number of Scientific Papers Published on Cancer By Month from 1964 to Present", ylab = "Total number of Papers", type = "l")
points (y = data.ts.new, x = time(data.ts.new), pch= as.vector(season(data.ts.new)))

# seasonal decomposition using loess
plot.stl <- function(..., xlab = "time") {
  mtext <- function(text, ...)
    graphics::mtext(if (text == "time") xlab else text, ...)
  plot.stl <- stats:::plot.stl
  environment(plot.stl) <- environment()
  plot.stl(...)
}

plot(stl ( data.ts.new, "periodic"))
swe = stl ( data.ts.new, "periodic")
plot.stl(swe, xlab = "Year")

#Using Lowess to look at the trend.
plot(lowess(diff(diff(sqrt(data.ts.new),lag = 6)), f  = 0.15))
plot(diff(diff(sqrt(data.ts.new)),lag = 6) - lowess(diff(diff(sqrt(data.ts.new)),lag = 6), f  = 0.1)$y)
points (data.ts.new - lowess(data.ts.new, f  = 0.1)$y, x = time(data.ts.new), pch= as.vector(season(data.ts.new)))


############################
#Transforming the data
############################

BoxCox.ar(data.ts.new,method = 'ols',lambda = seq(0.36,.65,0.01))
#suggests square root transformation


############################
#Ploting ACF and PACF of first and seasonal differences
##############################


# plot (diff(sqrt(data.ts.new)), main = "Plot of First Diff of Sq. Root")
# par(mfrow = c(1,2))
# acf (as.vector(diff(sqrt(data.ts.new))), main = "ACF of First Diff of Sq. Root")
# pacf (as.vector(diff(sqrt(data.ts.new))), main = "PACF of First Diff of Sq. Root")

# plot the first difference and acf, pacf
plot (diff(sqrt(data.ts.new)), ylab = "First Difference of Square Root Counts", xlab = "Year")
par(mfrow = c(1,2))
acf (as.vector(diff(sqrt(data.ts.new))), lag.max = 50, main = "")
pacf (as.vector(diff(sqrt(data.ts.new))),lag.max = 50, main = "")


#plot seasonal difference boxplots, then the first and seasonal difference ts
par(mfrow = c(1,3))

#boxplots
years = rep ( seq(1,12,1), 48)
boxplot ( diff(sqrt(data.ts.new)) ~ years[-1] , xlab = "Month", ylab = "First Difference of Square Root", main = "Before Seasonal Difference")
boxplot( diff(diff((sqrt(data.ts.new))),lag = 12) ~ years[-c(1:13)], xlab = "Month", ylab = "First and Seasonal Difference of Square Root", main = "After Seasonal Difference")

plot (diff(diff(sqrt(data.ts.new),lag = 12)), ylab = "First and Seasonal Difference of Square Root Counts", xlab = "Year", main = "First and Seasonal Difference")
abline( h = 0, col = "red")


#plot(diff(diff((data.ts.new)),lag = 12), type = 'l')

#plot the acf and pacf
par(mfrow = c(1,2))
acf(as.vector(diff(diff(sqrt(data.ts.new)),lag = 12)), lag.max = 40, main = "")
pacf(as.vector(diff(diff(sqrt(data.ts.new)),lag = 12)),lag.max = 60, main = "")




########################################
function to estimate all possible combinations of p, q, P, and Q
#####################################

get.aic = function ( coefs ) {
  p = coefs[1]
  q = coefs[2]
  P = coefs[3]
  Q = coefs[4]	
  mod = arima (diff(diff(sqrt(data.ts.new)),lag =6), order = c(p,0,q),seasonal=list(order=c(P,0,Q),period=6))
  #calcaicz = -2 * mod$loglik + 2*(p + q + P + Q + 2)
  return ( mod$aic )
  
}


aicvalues = vector ()
calcaic = vector ()
sarima.coeffs = matrix ( ncol = 5 )

for ( i in 1:5 ) {
  p = i-1
  for ( j in 1:5 ) {
    q = j-1
    for ( k in 1:4 ) {
      P = k-1
      for ( l in 1:4 ) {
        Q = l-1
        result <- try(get.aic (c(p,q,P,Q)));
        if(class(result) == "try-error") 
          next
        else{
          coefffs = cbind (rbind(c(p,q,P,Q)), rbind(result))
          sarima.coeffs  = rbind (sarima.coeffs, coefffs)
        }
      }
    }
  }
}


finalm7 = sarima.coeffs[-1,]

files = file("aic_values.txt")
write.table(finalm7,sep = ',',files)
close(files)

finalm[which ( finalm[,5] == min(finalm[,5])),]

###################
Picking and fitting the Model
###################


    p q P Q      AIC
198 2 2 2 3 1717.015
258 3 1 2 3 1717.473
364 4 3 3 1 1720.384
366 4 3 3 3 1721.304
361 4 3 2 2 1721.817
365 4 3 3 2 1721.954

model3 = arima (diff(diff(sqrt(data.ts.new),lag = 12)), order = c(2,0,2),seasonal=list(order=c(2,0,3),period=12))


Call:
arima(x = diff(diff(sqrt(data.ts.new), lag = 12)), order = c(2, 0, 2), seasonal = list(order = c(2, 
    0, 3), period = 12), io = c(24, 552))

0.5403  0.1825  -1.3925  0.4335  0.9938  -0.3069  -1.6927  0.8090  -0.0090

Coefficients:
         ar1     ar2      ma1     ma2    sar1     sar2     sma1    sma2     sma3 
      0.4854  0.2117  -1.3663  0.4006  0.7712  -0.9929  -1.5282  1.5405  -0.7428     
s.e.  0.0767  0.0620   0.0924  0.0893  0.0260   0.0014   0.0547  0.0431   0.0411    

sigma^2 estimated as 1.085:  log likelihood = -836.7,  aic = 1697.4


#############
Model Diagnostics
##############
resid = model3$residuals
par(mfrow=c(2,2))
spec(as.vector(resid), main = "Spectral Density of Residuals")
plot(as.vector(resid) ~ fitted(model3), main='Residuals Vs. Fitted', ylab = "Residual", xlab = "Fitted Value")
abline( h = 3, lty= 'dotted', col = 'blue')
abline (h = -3, lty= 'dotted', col = 'blue')
abline (h = 0, col = 'red')

acf(as.vector(resid), main='Sample ACF of Residuals', lag.max=20)
qqnorm(resid, main='Q-Q Plot of Estimated Residuals')
qqline(resid)

shapiro.test(resid)


Shapiro-Wilk normality test

data:  resid 
W = 0.9957, p-value = 0.1219





###########
Spectral Density Analysis
##########
k = kernel ( 'daniell',m = 8)
sp = spec ( diff(diff(sqrt(data.ts.new),lag = 12)), kernel = k, sub = '', ci.plot = T)


par(mfrow = c(2,1))
sp = spec (diff(diff(sqrt(data.ts.new)),lag = 12),span = 8, lty = 'dashed', ci.plot=T, col = "blue", main = "")
lines (sp$freq,sp$spec,lty='dashed')
f = seq(0.001,0.5,by = 0.001)
lines (f, ARMAspec(model = list(ma = c(-1.366258835 , 0.400618441 ), ar = c( 0.485441170 , 0.211741866), seasonal = list(sar=  c(0.771184427, -0.992886651), sma = c(-1.528197024 , 1.540478583, -0.742757029), period = 12)),freq = f,plot = F)$spec,lty="solid", col = "red")



