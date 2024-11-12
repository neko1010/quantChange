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

## plot the output 
plot(fit, ANOVA = T)
#plot(fit.NA, ANOVA = T)

######### Bayesian Change Point (for comparison) #############

##
#install.packages('bcp') 
library(bcp)

## include a 'year' variable
data.fill$year = as.numeric(format(as.Date(data.fill$system.time_start, format = "%b %e, %Y"), "%Y"))

## univariate
results.bcp = bcp(data.fill$mesic)
plot(results.bcp,main = "WRP Mesic Probabilities of Change - Yankee Fork",xlab = "Date", xaxlab = data.fill$year )

## with NAs  - DOESNT WORK
#results.bcp.na = bcp(data.NAs$mesic)
#plot(results.bcp,main = "WRP Mesic Probabilities of Change - Yankee Fork",xlab = "Date", xaxlab = data.NAs$date )


## NEXT STEPS - MULTIVARIATE - control for months as random effects, climate (drought) as fixed effect

#install.packages('Rbeast')
library(Rbeast)

results.beast = beast(data.fill$mesic, start = as.Date('2004-6-1'), deltat = "3 months") ##not quite right - missing 2016 data.... 
plot(results.beast, interactive = F)

## missing data
results.beast.na = beast(data.ts.NA) ## supposed to work but doesnt... https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf
results.beast.na = beast(data.ts.NA, start = c('2004, 1, 15'), deltat = 1/4, period = 365 )  
plot(results.beast.na, interactive = F)

## BSTS 
library(CausalImpact)

## data with drought index predictor
pdsi = read.csv('pdsi_Yankee.csv')

## remove 2016 to have corresponding values with the time series
pdsi$year = as.numeric(format(as.Date(pdsi$system.time_start, format = "%b %e, %Y"), "%Y"))
pdsi = subset(pdsi, year != 2016)

data.bsts=  cbind(data.fill$mesic, pdsi$pdsi)
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
