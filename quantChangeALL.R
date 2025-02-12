## install the package
#install.packages('bfast')

## load the required packages
library(bfast)

## read the data
data = read.csv('yankee.csv')

## BFAST can't handle missing values 
## Insert a row for July 2019
data[nrow(data) +1,] = c("July 1, 2019", NA)

## sort the df
data = data[order(as.Date(data$system.time_start, format = "%b %e, %Y")),]
rownames(data) = 1:nrow(data)

data$mesic = as.numeric(data$mesic)

## define a function to interpolate
replace_na_with_mean <- function(df) {
  for (i in 1:ncol(df)) {
    na_indices <- which(is.na(df[, i]))
    for (j in na_indices) {
      if (j == 1) {
        df[j, i] <- df[j + 1, i]
      } else if (j == nrow(df)) {
        df[j, i] <- df[j - 1, i]
      } else {
        df[j, i] <- mean(c(df[j - 1, i], df[j + 1, i]))
      }
    }
  }
  return(df)
}

## apply the interpolation function to fill NAs
data.fill = replace_na_with_mean(data)

## Embrace the NAs! insert some more for 2016
data.2016 = cbind(c("Jun 1, 2016", "July 1, 2016", "Aug 1, 2016", "Sep 1, 2016"),
                  c(rep(NA,4)))

colnames(data.2016) = c("system.time_start", "mesic")
data.NAs = rbind(data, data.2016)
data.NAs$date = as.Date(data.NAs$system.time_start, format = "%b %e, %Y")
data.NAs = data.NAs[order(data.NAs$date),]

## create a time series object
data.ts = ts(data.fill$mesic, frequency = 4)
data.ts.NA = ts(data.NAs$mesic, frequency = 4)

## assigning time units to observations
#tsp(data.ts.NA) = c(as.Date("2004-06-01"),as.Date("2020-09-01"), 4 )


plot(data.ts)
plot(data.ts.NA)

## apply the BFAST function
fit = bfast(data.ts, h = 0.15, season = "dummy")
#fit.na = bfast(data.ts.NA, h = 0.15, season = "dummy") ## can't handle missing data

## test the sensitivity to the prop of observations used to detect breakpoint; < 0.5
fit.05 = bfast(data.ts, h = 0.05, season = "dummy")
fit.25 = bfast(data.ts, h = 0.25, season = "dummy")
fit.35 = bfast(data.ts, h = 0.35, season = "dummy")
fit.45 = bfast(data.ts, h = 0.45, season = "dummy")

## plot the output 
plot(fit, ANOVA = T)

#plot(fit.05, ANOVA = T) ## no diffs
#plot(fit.25, ANOVA = T) ## no diffs
#plot(fit.35, ANOVA = T) ## no diffs
#plot(fit.45, ANOVA = T) ## no diffs
#plot(fit.NA, ANOVA = T) ## doesn't work

######### Bayesian Change Point (for comparison) #############

##
#install.packages('bcp') 
library(bcp)

## include a 'year' variable
data.fill$year = as.numeric(format(as.Date(data.fill$system.time_start, format = "%b %e, %Y"), "%Y"))

## univariate
results.bcp = bcp(data.fill$mesic)
plot(results.bcp,main = "WRP Mesic Probabilities of Change - Yankee Fork",xlab = "Date", xaxlab = data.fill$year )

summary(results.bcp)

## Multivariate - use PDSI as a predictor with time

## data with drought index predictor
pdsi = read.csv('pdsi_Yankee.csv')

## remove 2016 to have corresponding values with the time series
pdsi$year = as.numeric(format(as.Date(pdsi$system.time_start, format = "%b %e, %Y"), "%Y"))
pdsi.no2016 = subset(pdsi, year != 2016)

data.multiv = cbind(data.fill, pdsi.no2016$pdsi)
results.bcp.multiv = bcp(data.multiv$mesic, data.multiv$`pdsi.no2016$pdsi`)
plot(results.bcp.multiv,main = "WRP Mesic Probabilities of Change - Yankee Fork Multivariate",xlab = "Date", xaxlab = data.fill$year )

summary(results.bcp.multiv)
 
## with NAs  - DOESNT WORK
#results.bcp.na = bcp(data.NAs$mesic)
#plot(results.bcp,main = "WRP Mesic Probabilities of Change - Yankee Fork",xlab = "Date", xaxlab = data.NAs$date )

## Test sensitivity of the w0 and p0 params - quant diffs minimal, no qualitative diffs
#results.bcp.w1 = bcp(data.fill$mesic, w0 = 0.1)
#results.bcp.w3 = bcp(data.fill$mesic, w0 = 0.3)
#results.bcp.w4 = bcp(data.fill$mesic, w0 = 0.4)
#results.bcp.w5 = bcp(data.fill$mesic, w0 = 0.5)
#results.bcp.w6 = bcp(data.fill$mesic, w0 = 0.6)
#results.bcp.w7 = bcp(data.fill$mesic, w0 = 0.7)
#results.bcp.w8 = bcp(data.fill$mesic, w0 = 0.8)
#results.bcp.w9 = bcp(data.fill$mesic, w0 = 0.9)
#
#results.bcp.p1 = bcp(data.fill$mesic, p0 = 0.1)
#results.bcp.p3 = bcp(data.fill$mesic, p0 = 0.3)
#results.bcp.p4 = bcp(data.fill$mesic, p0 = 0.4)
#results.bcp.p5 = bcp(data.fill$mesic, p0 = 0.5)
#results.bcp.p6 = bcp(data.fill$mesic, p0 = 0.6)
#results.bcp.p7 = bcp(data.fill$mesic, p0 = 0.7)
#results.bcp.p8 = bcp(data.fill$mesic, p0 = 0.8)
#results.bcp.p9 = bcp(data.fill$mesic, p0 = 0.9)
#
#plot(results.bcp.w1, main = "WRP Mesic Probabilities of Change - w0 = 0.1",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w3, main = "WRP Mesic Probabilities of Change - w0 = 0.3",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w4, main = "WRP Mesic Probabilities of Change - w0 = 0.4",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w5, main = "WRP Mesic Probabilities of Change - w0 = 0.5",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w6, main = "WRP Mesic Probabilities of Change - w0 = 0.6",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w7, main = "WRP Mesic Probabilities of Change - w0 = 0.7",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w8, main = "WRP Mesic Probabilities of Change - w0 = 0.8",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.w9, main = "WRP Mesic Probabilities of Change - w0 = 0.9",xlab = "Date", xaxlab = data.fill$year )
#
#plot(results.bcp.p1, main = "WRP Mesic Probabilities of Change - p0 = 0.1",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p3, main = "WRP Mesic Probabilities of Change - p0 = 0.3",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p4, main = "WRP Mesic Probabilities of Change - p0 = 0.4",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p5, main = "WRP Mesic Probabilities of Change - p0 = 0.5",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p6, main = "WRP Mesic Probabilities of Change - p0 = 0.6",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p7, main = "WRP Mesic Probabilities of Change - p0 = 0.7",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p8, main = "WRP Mesic Probabilities of Change - p0 = 0.8",xlab = "Date", xaxlab = data.fill$year )
#plot(results.bcp.p9, main = "WRP Mesic Probabilities of Change - p0 = 0.9",xlab = "Date", xaxlab = data.fill$year )

#install.packages('Rbeast')
library(Rbeast)

results.beast = beast(data.fill$mesic, start = as.Date('2004-6-1'), deltat = "3 months", dump.ci = T) ##not quite right - missing 2016 data.... 
plot(results.beast, interactive = F)

## missing data
results.beast.na = beast(as.numeric(data.ts.NA), start = as.Date('2004-1-1'), deltat = "3 months") 
plot(results.beast.na, interactive = F) ## missing data leads to higher likelihood of negative slopes, and an extra breakpoint estimated in the trend

## Parameters 
results.beast.na.unifPrec = beast(as.numeric(data.ts.NA), start = as.Date('2004-1-1'), deltat = "3 months", precPriorType = 'uniform') 
plot(results.beast.na.unifPrec, interactive = F) 

results.beast.na.compPrec = beast(as.numeric(data.ts.NA), start = as.Date('2004-1-1'), deltat = "3 months", precPriorType = 'componentwise') 
plot(results.beast.na.compPrec, interactive = F) 

## BSTS 
library(CausalImpact)

## Include PDSI as a predictor
data.bsts=  cbind(data.fill$mesic, pdsi.no2016$pdsi)
data.bsts.time =  zoo(cbind(data.fill$mesic, pdsi$pdsi), as.Date(data.fill$system.time_start,format = "%b %e, %Y")) ## includes missing data period

## define pre and post restoration periods
pre.period = c(1,32)
post.period = c(33, 64)

## with date
pre.period.date = as.Date(c("Jun 1, 2004", "Sep 1,2011"), format = "%b %e, %Y")
post.period.date = as.Date(c("Jun 1, 2012", "Sep 1,2020"), format = "%b %e, %Y")

impact = CausalImpact(data.bsts, pre.period, post.period, model.args = list(nseasons = 4))
impact.time = CausalImpact(data.bsts.time, pre.period.date, post.period.date, model.args = list(nseasons = 4))
plot(impact)
plot(impact.time)

summary(impact)
summary(impact.time) ## nearly identical results AND includes period of missing data...


library(terra)
ras = rast("Yankee.tif")
plot(ras)
