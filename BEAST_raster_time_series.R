###################################
# Applying BEAST to a raster time series
# Author: Carolyn Koehn
# Last modified: 4/8/25
###################################

# load packages
library(terra)
library(sf)
library(Rbeast)

# load data (June-Sep 2004-2020 with data gaps)
ras <- rast("Yankee.tif")

# get shapefile of area of interest
aoi <- st_read("yankee.shp") %>%
  st_transform(crs(ras))

############################################
# Data Preparation

# add time stamps to raster
names = names(ras)
dates = sapply(names, function(name){
  as.Date(paste0(unlist(strsplit(name, "_"))[4],"01"),"%Y%m%d")
})
dates.sort = order(dates)
ras.sort = ras[[dates.sort]]

# assign times to layers (used in BEAST)
time(ras.sort) <- as.Date(dates[dates.sort])

# convert to array
array.input <- as.array(ras.sort)

##############################################
# apply BEAST to 3D time series
# "pseudo-BFAST" method
beast.output <- beast123(Y = array.input,
                         metadata = 
                           list(whichDimIsTime = 3, # the layer dimension (3rd dimension) corresponds to time
                                startTime = as.Date('2004-1-1'),
                                deltaTime = "3 months", # time between data points
                                period = "1 year"), # length of seasonal period (we expect a yearly cycle)
                         season = 'harmonic',
                         mcmc = list(seed = 101), # set seed for +/- reproducible results
                         extra = list(dumpInputData = TRUE)) # return the data used by BEAST (the regular time series generated from our irregular time series) to illustrate what is occurring under the hood

# irregular time series method
beast.output.irreg <- beast123(Y = array.input,
                         metadata = 
                           list(whichDimIsTime = 3, # the layer dimension (3rd dimension) corresponds to time
                                time = time(ras.sort), # since time series is irregular, provide times here
                                startTime = time(ras.sort)[1],
                                deltaTime = "1 month", # time between data points
                                period = "1 year", # length of seasonal period (we expect a yearly cycle) 
                                sorder.minmax = 2), # we expect only one annual max/min
                         season = 'harmonic',
                         mcmc = list(seed = 101), # set seed for +/- reproducible results
                         extra = list(dumpInputData = TRUE)) # return the data used by BEAST (the regular time series generated from our irregular time series) to illustrate what is occurring under the hood

#############################################
# Reading the output

## Our data had 11 rows and 26 columns in the raster, with 63 time points
## BEAST filled in our irregular time series, adding NAs for missing months
## Here is an example of the time series generated for the top left pixel:
print(beast.output$data[1,1,])

## see summary plot for one pixel
plot(beast.output, index = c(15, 7))

# Measures of model fit by pixel
# R2
par(mfrow=c(1,2))
plot(rast(beast.output$R2,
          crs=crs(ras.sort),
          extent=ext(ras.sort)), # convert matrix to raster
     main = "R2 of model fit")
# RMSE
plot(rast(beast.output$RMSE,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "RMSE of model fit")

# Time series descriptions
## Number of change points (median over all models)
par(mfrow = c(1,2))
plot(rast(beast.output$trend$ncp_median,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Median number of\ntrend change points")
plot(rast(beast.output$season$ncp_median,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Median number of\nseasonal change points")

## Probability of the median number of change points
get_Pr_for_ncp <- function(ncpPr, ncp) {
  return(ncp.Pr[[1]][ncp+1])
}

ncpPr_trend_median <- matrix(data = mapply(FUN = get_Pr_for_ncp,
                                           apply(beast.output$trend$ncpPr, MARGIN=c(1,2), list),
                                           beast.output$trend$ncp_median),
                             nrow = nrow(beast.output$trend$ncp_median),
                             ncol = ncol(beast.output$trend$ncp_median),
                             byrow = FALSE)

ncpPr_season_median <- matrix(data = mapply(FUN = get_Pr_for_ncp,
                                           apply(beast.output$season$ncpPr, MARGIN=c(1,2), list),
                                           beast.output$season$ncp_median),
                             nrow = nrow(beast.output$season$ncp_median),
                             ncol = ncol(beast.output$season$ncp_median),
                             byrow = FALSE)
  
par(mfrow = c(1,2))
plot(rast(ncpPr_trend_median,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Probability of\nmedian number of\ntrend change points")
plot(rast(ncpPr_season_median,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Probability of\nmedian number of\nseasonal change points")

## Largest magnitude of change in trend
get_largest_magn_change_stats <- function(cpAbruptChange, ncp, cpPr, cp) {
  
}

cpMagn_trend <- matrix(data=NA, 
                       nrow=nrow(beast.output$trend$cpAbruptChange), 
                       ncol=ncol(beast.output$trend$cpAbruptChange))
for(i in 1:nrow(beast.output$trend$cpAbruptChange)) {
  for(j in 1:ncol(beast.output$trend$cpAbruptChange)) {
    if(beast.output$trend$ncp_median[i, j] > 0) {
      highest_magn_trend <- which.max(abs(beast.output$trend$cpAbruptChange[i, j, ]))
      cpMagn_trend[i,j] <- beast.output$trend$cpAbruptChange[i, j, highest_magn_trend]
    } else {
      cpMagn_trend[i,j] <- NA
    }
  }
}

### Probability of largest change
cpMagnProb_trend <- matrix(data=NA, 
                           nrow=nrow(beast.output$trend$cpPr), 
                           ncol=ncol(beast.output$trend$cpPr))
for(i in 1:nrow(beast.output$trend$cpPr)) {
  for(j in 1:ncol(beast.output$trend$cpPr)) {
    if(beast.output$trend$ncp_median[i, j] > 0) {
      highest_magn_trend <- which.max(abs(beast.output$trend$cpAbruptChange[i, j, ]))
      cpMagnProb_trend[i,j] <- beast.output$trend$cpPr[i, j, highest_magn_trend]
    } else {
      cpMagnProb_trend[i,j] <- NA
    }
  }
}

### Time of largest magnitude of change in trend
cpMagnTime_trend <- matrix(data=NA, 
                           nrow=nrow(beast.output$trend$cp), 
                           ncol=ncol(beast.output$trend$cp))
for(i in 1:nrow(beast.output$trend$cp)) {
  for(j in 1:ncol(beast.output$trend$cp)) {
    if(beast.output$trend$ncp_median[i, j] > 0) {
      highest_magn_trend <- which.max(abs(beast.output$trend$cpAbruptChange[i, j, ]))
      cpMagnTime_trend[i,j] <- beast.output$trend$cp[i, j, highest_magn_trend]
    } else {
      cpMagnProb_trend[i,j] <- NA
    }
  }
}

### Largest magnitude in trend plots
par(mfrow = c(1,3))
plot(rast(cpMagn_trend,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Panel A\nLargest magnitude\nof change in trend",
     mar = c(2.6, 4.6, 3.6, 6.6),
     col = viridisLite::viridis(50))
plot(st_geometry(aoi), border="black", lwd = 3, add=T)
plot(rast(cpMagnProb_trend,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Panel B\nProbability of largest\nchange magnitude in trend",
     mar = c(2.6, 4.6, 3.6, 6.6),
     col = viridisLite::viridis(50))
plot(st_geometry(aoi), border="black", lwd = 3, add=T)
plot(rast(cpMagnTime_trend,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Panel C\nTime of largest change\nmagnitude in trend",
     mar = c(2.6, 4.6, 3.6, 6.6),
     col = viridisLite::viridis(50))
plot(st_geometry(aoi), border="black", lwd = 3, add=T)

## Trend component from all pixels
par(mfrow=c(1,1))
plot(x = beast.output$time,
     xlab = "Date",
     y = beast.output$trend$Y[1,1,],
     ylab = "Time Series Trend",
     type="l", main = "Trend Component for all Pixels",
     ylim = c(0,0.8), col="gray80") # y-limits should be adjusted for your results
for(i in 1:nrow(beast.output$trend$Y)) {
  for(j in 1:ncol(beast.output$trend$Y)) {
    lines(x = beast.output$time,
          y = beast.output$trend$Y[i,j,], 
          col="gray80")
  }
}




## Largest magnitude of change in seasonal
cpMagn_season <- matrix(data=NA, 
                       nrow=nrow(beast.output$season$cpAbruptChange), 
                       ncol=ncol(beast.output$season$cpAbruptChange))
for(i in 1:nrow(beast.output$season$cpAbruptChange)) {
  for(j in 1:ncol(beast.output$season$cpAbruptChange)) {
    if(beast.output$season$ncp_median[i, j] > 0) {
      highest_magn_season <- which.max(abs(beast.output$season$cpAbruptChange[i, j, ]))
      cpMagn_season[i,j] <- beast.output$season$cpAbruptChange[i, j, highest_magn_season]
    } else {
      cpMagn_season[i,j] <- NA
    }
  }
}

### Probability of largest change
cpMagnProb_season <- matrix(data=NA, 
                           nrow=nrow(beast.output$season$cpPr), 
                           ncol=ncol(beast.output$season$cpPr))
for(i in 1:nrow(beast.output$season$cpPr)) {
  for(j in 1:ncol(beast.output$season$cpPr)) {
    if(beast.output$season$ncp_median[i, j] > 0) {
      highest_magn_season <- which.max(abs(beast.output$season$cpAbruptChange[i, j, ]))
      cpMagnProb_season[i,j] <- beast.output$season$cpPr[i, j, highest_magn_season]
    } else {
      cpMagnProb_season[i,j] <- NA
    }
  }
}

### Time of largest magnitude of change in season
cpMagnTime_season <- matrix(data=NA, 
                           nrow=nrow(beast.output$season$cp), 
                           ncol=ncol(beast.output$season$cp))
for(i in 1:nrow(beast.output$season$cp)) {
  for(j in 1:ncol(beast.output$season$cp)) {
    if(beast.output$season$ncp_median[i, j] > 0) {
      highest_magn_season <- which.max(abs(beast.output$season$cpAbruptChange[i, j, ]))
      cpMagnTime_season[i,j] <- beast.output$season$cp[i, j, highest_magn_season]
    } else {
      cpMagnProb_season[i,j] <- NA
    }
  }
}

### Largest magnitude in season plots
par(mfrow = c(1,3))
plot(rast(cpMagn_season,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Largest magnitude\nof change in season")
plot(rast(cpMagnProb_season,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Probability of\nlargest magnitude\nof change in season")
plot(rast(cpMagnTime_season,
          crs=crs(ras.sort),
          extent=ext(ras.sort)),
     main = "Time of\nlargest magnitude\nof change in season")

## Seasonal component from all pixels
par(mfrow=c(1,1))
plot(x = beast.output$time,
     xlab = "Date",
     y = beast.output$season$Y[1,1,],
     ylab = "Time Series Seasonal",
     type="l", main = "Seasonal Component for all Pixels",
     ylim = c(-0.6,0.6), col="gray80") # y-limits should be adjusted for your results
for(i in 1:nrow(beast.output$season$Y)) {
  for(j in 1:ncol(beast.output$season$Y)) {
    lines(x = beast.output$time,
          y = beast.output$season$Y[i,j,], 
          col="gray80")
  }
}
