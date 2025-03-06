###################################
# Applying BEAST to a single time series
# Author: Carolyn Koehn
# Last modified: 3/6/25
###################################

# load packages
library(Rbeast)

# read the data
data <- read.csv('yankee.csv')

# apply datetime format to data
data$system.time_start <- as.Date(data$system.time_start, format = "%b %e, %Y")

# apply BEAST
# using beast.irreg because we only have summer data in this example
# if your data is evenly spaced in time, use beast function
beast.output <- beast.irreg(y = data$mesic, # vector for time series
                            time = data$system.time_start, # vector of times
                            deltat = "1 month", # time between observations
                            season = "harmonic", # we expect an annual peak and low
                            period = "1 year",  # we expect an annual peak and low
                            sorder.minmax = 1, # we expect one annual peak and low
                            mcmc.seed = 101) # set seed for +/- reproducible results

# plot summary results
plot(beast.output)
