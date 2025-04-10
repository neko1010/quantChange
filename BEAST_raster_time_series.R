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

## Descriptive stats for largest magnitude of change in trend
get_largest_magn_change_stats <- function(cpAbruptChange, ncp, cpPr, cp) {
  stats.to.return <- c(cpMagn = NA,
                       cpMagnPr = NA,
                       cpMagnTime = NA)
  if(ncp > 0) {
    highest.magn.trend <- which.max(abs(cpAbruptChange[[1]]))
    
    stats.to.return["cpMagn"] <- cpAbruptChange[[1]][highest.magn.trend]
    stats.to.return["cpMagnPr"] <- cpPr[[1]][highest.magn.trend]
    stats.to.return["cpMagnTime"] <- cp[[1]][highest.magn.trend]
  }
  return(stats.to.return)
}

largest_magn_trend_change <- mapply(FUN = get_largest_magn_change_stats,
                 apply(beast.output$trend$cpAbruptChange, MARGIN = c(1,2), list),
                 beast.output$trend$ncp_median,
                 apply(beast.output$trend$cpPr, MARGIN = c(1,2), list),
                 apply(beast.output$trend$cp, MARGIN = c(1,2), list))

### Magnitude of largest trend change
cpMagn_trend <- matrix(data = largest_magn_trend_change["cpMagn",],
                       nrow = nrow(beast.output$trend$cpAbruptChange),
                       ncol = ncol(beast.output$trend$cpAbruptChange),
                       byrow = FALSE)

### Probability of largest trend change
cpMagnProb_trend <- matrix(data = largest_magn_trend_change["cpMagnPr",],
                       nrow = nrow(beast.output$trend$cpAbruptChange),
                       ncol = ncol(beast.output$trend$cpAbruptChange),
                       byrow = FALSE)

### Time of largest trend change
cpMagnTime_trend <- matrix(data = largest_magn_trend_change["cpMagnTime",],
                       nrow = nrow(beast.output$trend$cpAbruptChange),
                       ncol = ncol(beast.output$trend$cpAbruptChange),
                       byrow = FALSE)

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
apply(beast.output$trend$Y,
      MARGIN = c(1,2),
      function(y) {
        lines(x = beast.output$time,
              y = y, 
              col="gray80")
        return()
      })


## Largest seasonal change magnitude
largest_magn_season_change <- mapply(FUN = get_largest_magn_change_stats,
                                    apply(beast.output$season$cpAbruptChange, MARGIN = c(1,2), list),
                                    beast.output$season$ncp_median,
                                    apply(beast.output$season$cpPr, MARGIN = c(1,2), list),
                                    apply(beast.output$season$cp, MARGIN = c(1,2), list))

### Magnitude of largest seasonal change
cpMagn_season <- matrix(data = largest_magn_season_change["cpMagn",],
                       nrow = nrow(beast.output$season$cpAbruptChange),
                       ncol = ncol(beast.output$season$cpAbruptChange),
                       byrow = FALSE)

### Probability of largest seasonal change
cpMagnProb_season <- matrix(data = largest_magn_season_change["cpMagnPr",],
                           nrow = nrow(beast.output$season$cpAbruptChange),
                           ncol = ncol(beast.output$season$cpAbruptChange),
                           byrow = FALSE)

### Time of largest seasonal change
cpMagnTime_season <- matrix(data = largest_magn_season_change["cpMagnTime",],
                           nrow = nrow(beast.output$season$cpAbruptChange),
                           ncol = ncol(beast.output$season$cpAbruptChange),
                           byrow = FALSE)

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
apply(beast.output$season$Y,
      MARGIN = c(1,2),
      function(y) {
        lines(x = beast.output$time,
              y = y, 
              col="gray80")
        return()
      })
