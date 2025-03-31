###################################
# Applying BSTS to a raster time series
# Author: Juan C Rojas L
# Last modified: 02/22/25
###################################
#The need packages
library(bfast)
library(terra)
library(gtools)
library(CausalImpact)
library(terra)
library(zoo)

setwd('C:\\Users\\JUANCAMILOROJASL\\Desktop\\BSTS\\quantChange-main\\') #set the location
ras = rast('Yankee.tif')
plot(ras)
dim(ras)
## insert an entire month of NAs for July 2019
july19 = rast(ras[[1]], vals = NA, names = c('1_1_040029_201907_mesic'))
## add this empty rast to the full stack
ras = c(ras, july19)
print(ras)
## add time stamps to each
names = names(ras)
dates = sapply(names, function(name){
  as.Date(paste0(unlist(strsplit(name, "_"))[4],"01"),"%Y%m%d")
})
## sort by dates
dates.sort = order(dates)
ras.sort = ras[[dates.sort]]
time(ras.sort) 
dim(ras.sort)
##################################################################### 
############################################################
##################################################
### Uploading data.
data <- as.data.frame(ras.sort)
pdsi = read.csv('pdsi_Yankee.csv')
## remove 2016 to have corresponding values with the time series
pdsi$year = as.numeric(format(as.Date(pdsi$system.time_start, format = "%b %e, %Y"), "%Y"))
pdsi.no2016 = subset(pdsi, year != 2016)
pdsi.repeated<-(pdsi.no2016$pdsi)
pdsi_matrix <- matrix(rep(pdsi.repeated, each = 286), ncol = length(pdsi.repeated))
dim(pdsi_matrix)
dim(data)
#I've created the matrix since I need a pdsi value per pixel over time. 
#Considering that PDSI has a value per year over that area, what I made was repeat the values in a matrix. 
#########################
# Initialize an empty list to store combined zoo objects
combined_zoo <- list()
all_values <- c()# Initialize an empty vector to store values
data[] <- lapply(seq_along(data), function(i) {
  x <- data[[i]]  # Extract the current column
  # If the entire column is NA, replace it with the mean of the previous and next columns
  if (all(is.na(x))) {
    if (i > 1 & i < ncol(data)) {  
      x[] <- rowMeans(cbind(data[[i - 1]], data[[i + 1]]), na.rm = TRUE)
      print(x[])
    } else if (i == 1) {  # If it's the first column, use only the next column
      x[] <- data[[i + 1]]
    } else if (i == ncol(data)) {  # If it's the last column, use only the previous column
      x[] <- data[[i - 1]]
    }
  } else {
    x[is.na(x)] <- mean(x, na.rm = TRUE)  # Replace individual NAs with the column mean
  }
  
  return(x)
})

all_values <- matrix(NA, nrow = nrow(data), ncol = 64)  
max_values_ci <- matrix(NA, nrow = nrow(data), ncol = 64) 
# Loop through each row of the data
for (i in 1:nrow(data)) {
  column_mesic <- as.numeric(data[i,])
  column_pdsi<- pdsi_matrix[i,]
  combined_zoo<- zoo(cbind(column_mesic, column_pdsi))
  pre.period <- c(1, 32)
  post.period <- c(33, 64)
  impact <- CausalImpact(combined_zoo, pre.period, post.period, model.args = list(nseasons = 4))
  all_values[i,] <- impact$series$point.effect# Return the point effect of the series
  for(j in 1:ncol(all_values)) {
    max_values_ci[i,j] <- max(all_values[i,], na.rm = TRUE)    # Get the highest value in column 'i'
  }
}
dim(all_values)
magn <- max_values_ci
point_effect_matrix <- matrix(magn, nrow = nrow(ras.sort), ncol = ncol(ras.sort), byrow = TRUE)
point_effect_raster <- rast(point_effect_matrix)
ext(point_effect_raster) <- ext(ras.sort)
crs(point_effect_raster) <- crs(ras.sort) 
plot(point_effect_raster, main = "Largest Magnitude of Change in Trend", range = c(0, 0.6))
print(impact)
######################### Custom model 
all_values_fl <- matrix(NA, nrow = nrow(data), ncol = 64)  
max_values <- matrix(NA, nrow = nrow(data), ncol = 64) 
for (i in 1:nrow(data)) {
  column_mesic <- as.numeric(data[i,])
  column_pdsi<- pdsi_matrix[i,]
  post.period <- c(33, 64)
  post.period.response <- column_mesic[post.period[1] : post.period[2]]
  column_mesic[post.period[1] : post.period[2]] <- NA
  ss <- AddLocalLevel(list(), column_mesic)
  ss <- AddSeasonal(ss, column_mesic, nseasons = 4)
  bsts.model <- bsts(column_mesic ~ column_pdsi, ss, niter = 1000)
  impact <- CausalImpact(bsts.model = bsts.model,
                         post.period.response = post.period.response)
  all_values_fl[i,] <- impact$series$point.effect# Return the point effect of the series
  for(j in 1:ncol(all_values_fl)) {
    max_values[i,j] <- max(all_values_fl[i,], na.rm = TRUE)    # Get the highest value in column 'i'
  }
}
summary(impact)
#print(max_values)
magn <- max_values
point_effect_matrix <- matrix(magn, nrow = nrow(ras.sort), ncol = ncol(ras.sort), byrow = TRUE)
point_effect_raster <- rast(point_effect_matrix)
ext(point_effect_raster) <- ext(ras.sort)
crs(point_effect_raster) <- crs(ras.sort) 
# Adjust margins to allow space for long titles
par(mfrow = c(1, 3), mar = c(6, 4, 8, 2))  #
# Plot the largest magnitude of change in trend
plot(ras.sort[[1]], 
     main = "Panel A\nMesic Vegetation Proportion Sep 2004",
     zlim = c(0, 0.6),  cex.main = 1)
plot(ras.sort[[64]], 
     main = "Panel B\nMesic Vegetation Proportion Sep 2020",
     zlim = c(0, 0.6),  cex.main = 1)
plot(point_effect_raster, 
     main = "Panel C\nLargest Magnitude of Change in Trend", 
     zlim = c(0, 0.6),
     cex.main = 1)
# Reset layout to default
par(mfrow = c(1, 1))  