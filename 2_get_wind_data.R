# This script extracts the wind speed and direction for every survey day in the study area
# The data are at a 0.5 degree resolution

library(rWind)
# this creates an individual file for each day
wind_vals <- wind.dl_2(time = surv_dates, lon1 = min(nnpp$longitude), lon2 = max(nnpp$longitude), 
                       lat1 = min(nnpp$latitude), lat2 = max(nnpp$latitude), type = 'csv')

# we can see the list of files created above
wind_files <- list.files( "C:/Users/ebrowning/Dropbox/PhD/BCT_internship/analyses/data/wind_data", pattern = 'wind_*')

# now we need to bind all of those files together into one csv
setwd("C:/Users/ebrowning/Dropbox/PhD/BCT_internship/analyses/data/wind_data")
# read in the csv files created and bind to one dataframe
# also make a stack of speed and direction
wind_data <- data.frame()
wind_stack <- stack()
for (i in 1:length(wind_files)){
  file_i = read.csv(wind_files[i])
  file_name = strsplit(wind_files[i], split = "[.]")[[1]][1]
  #wind_data <- rbind(wind_data, file_i)
  wind_ras_i <- wind2raster(file_i)
  names(wind_ras_i)[[1]] <- paste(file_name, names(wind_ras_i)[[1]], sep = "_")
  names(wind_ras_i)[[2]] <- paste(file_name, names(wind_ras_i)[[2]], sep = "_")
  wind_stack <- stack(wind_stack, wind_ras_i)
}

write.csv(wind_data, 'all_wind_data.csv')
