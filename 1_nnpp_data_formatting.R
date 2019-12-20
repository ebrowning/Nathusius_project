#------------------------------#
# This script does some formatting of the dataset produced by the Access query 
 ## At the bottom of this script are some visulisations and explorations of the data 
#------------------------------#

setwd("~/Dropbox/PhD/BCT_internship/analyses/data") # change this to yours 
#------------------------------
# Data formating
#------------------------------
nnpp <- read.csv('nnpp_2014_18.csv', na.strings = c("", "NA")) # fill in blanks with NAs
head(nnpp)


# set -99 temp values to NA
nnpp$StartTemp_2 <- nnpp$StartTemp
nnpp$StartTemp_2[nnpp$StartTemp_2 < 0] <- NA

nnpp$EndTemp_2 <- nnpp$EndTemp
nnpp$EndTemp_2[nnpp$EndTemp_2 < 0] <- NA

# create columns for start and end dates and 'bat time' in unix time
#check 'strptime' help page for codes for dates
nnpp$surveystart_date_time <- paste(nnpp$SurveyDateStart, nnpp$SurveyTimeStart, sep = " ")
nnpp$surveyend_date_time <- paste(nnpp$SurveyDateEnd, nnpp$SurveyTimeEnd, sep = " ")
nnpp$surveytime_start_unix <- as.numeric(as.POSIXct(nnpp$surveystart_date_time, format="%d-%b-%y %H:%M:%S")) 
nnpp$surveytime_end_unix <- as.numeric(as.POSIXct(nnpp$surveyend_date_time, format="%d-%b-%y %H:%M:%S")) 
nnpp$bat_date_time <- paste(nnpp$BatDate, nnpp$BatTime, sep = " ")
nnpp$battime_unix <- as.numeric(as.POSIXct(nnpp$bat_date_time, format="%d-%b-%y %H:%M:%S")) 

#calculate survey length in minutes
nnpp$surveylength <- (nnpp$surveytime_end_unix - nnpp$surveytime_start_unix)/60

# Get land cover data
## 
library(raster)
## extract land cover data - dominant target class from LCM data
lcm <- raster('lcm2015_gb_1km_dominant_target_class.tif')

coord <- cbind(nnpp$Easting_1km, nnpp$Northing_1km)
lcm_vals <- extract(lcm, coord)
#group classess into aggregates, based mostly on LCM docmentation, but also have grouped saltwater and coastal together as coastal
nnpp$lcm_class <- lcm_vals
nnpp$lcm_class_names <- ifelse(nnpp$lcm_class == 1 | nnpp$lcm_class == 2, 'woodland',
                               ifelse(nnpp$lcm_class == 3, 'arable', 
                                      ifelse(nnpp$lcm_class == 4 | nnpp$lcm_class == 5, 'grassland', 
                                             ifelse(nnpp$lcm_class == 8, 'fen_marsh', 
                                                    ifelse(nnpp$lcm_class == 9 , 'heather',
                                                           ifelse(nnpp$lcm_class == 14, 'freshwater',
                                                                  ifelse(nnpp$lcm_class == 13 | nnpp$lcm_class == 16 | nnpp$lcm_class == 18 | nnpp$lcm_class == 19, 'coastal',
                                                                         ifelse(nnpp$lcm_class == 20 | nnpp$lcm_class == 21, 'urban', 'unknown'))))))))

# Calculate time since sunset that a bat was caught 
## Need to transform eastings/northings to long/lat
library(maptools)
library(rgdal)
### shortcuts
ukgrid <- "+init=epsg:27700"
latlong <- "+init=epsg:4326"
### Create coordinates variable
coords <- cbind(Easting = nnpp$Easting,
                Northing = nnpp$Northing)
### Create the SpatialPointsDataFrame
dat_SP <- SpatialPointsDataFrame(coords,
                                 data = nnpp,
                                 proj4string = CRS(ukgrid))
### Convert
dat_SP_LL <- spTransform(dat_SP, CRS(latlong))
# rename relevant columns
dat_SP_LL@data$longitude <- coordinates(dat_SP_LL)[, 1]
dat_SP_LL@data$latitude <- coordinates(dat_SP_LL)[, 2]
# bind the long/lat columns to main nnpp df
nnpp2 <- cbind(nnpp, longitude = dat_SP_LL@data$longitude, latitude = dat_SP_LL@data$latitude)

# remove rows where surveys weren't actually carried out 
nnpp3 <- nnpp2[!is.na(nnpp2$SurveyDateStart), ]

# now calculate the sunset time for each day
surveydate = as.POSIXct(nnpp3$SurveyDateStart, format="%d-%b-%y")
locs = cbind(nnpp3$longitude, nnpp3$latitude)
sunset = crepuscule(locs, surveydate, direction = 'dusk', solarDep = 1, POSIXct.out = T)
nnpp3$sunset = sunset$time
#convert to unix
nnpp3$sunset_unix = strftime(nnpp3$sunset, '%s')
# calculate difference in time bat was caught and sunset time (time in minutes)
nnpp3$battime_since_sunset = as.integer((nnpp3$battime_unix - as.numeric(nnpp3$sunset_unix))/60)

# also calculate week of the year each survey was carried out using the survey start date
nnpp3$week_of_year <-as.numeric(strftime(as.POSIXct(nnpp3$SurveyDateStart, format="%d-%b-%y"), format = '%V'))
# year
nnpp3$year <-as.numeric(strftime(as.POSIXct(nnpp3$SurveyDateStart, format="%d-%b-%y"), format = '%Y'))
#month
nnpp3$month <-as.numeric(strftime(as.POSIXct(nnpp3$SurveyDateStart, format="%d-%b-%y"), format = '%m'))

head(nnpp3)
# remove years 2011-2013 as very few recrds
nnpp4 <- nnpp3[nnpp3$year > 2013, ]

# round up coordinates to the nearest 1000 as lowest common denominator is 1km resolution
nnpp4$Easting_1km <- round(nnpp4$Easting, -2)
nnpp4$Northing_1km <- round(nnpp4$Northing, -2)

write.csv(nnpp4, 'nnpp_2014_18.csv')

#----------------------------- End of data formting ----------------------------------- #
# ------------------------------
# ------------------------------
# --------------------- Some data exploration 

## This section is a bit messy

nnpp <- read.csv('nnpp_2014_18.csv', na.strings = c("", "NA")) # fill in blanks with NAs
head(nnpp)

# how many surveys carried out
length(unique(nnpp3$SurveyID))

# how many records per year
table(nnpp3$year)

library(ggplot2)
ggplot(nnpp, aes(x = lcm_class_names))+ geom_histogram(stat = 'count')

ggplot(nnpp, aes(x = month, fill = Sex)) +geom_histogram(stat = 'count')
ggplot( )
nats <- nnpp[nnpp$SpeciesID == 23, ]
ggplot(nats, aes(x = month, fill = Age)) + geom_histogram(stat = 'count')


#---- visualise the data ----#
library(ggplot2)
ggplot(nnpp4, aes(month, BatNumber)) + geom_boxplot()
ggplot(nnpp4, aes(month)) +geom_histogram(stat = 'count')

ggplot(nnpp4, aes(Easting, Northing, size=BatNumber, col=week_of_year)) + geom_point() +facet_wrap(LatinName~year)

ggplot(nnpp4, aes(Easting, Northing, size=BatNumber)) + geom_point() +facet_wrap(~year)

ggplot(nnpp4, aes(TrapEasting, TrapNorthing, col=LatinName)) + geom_point() +facet_wrap(~year)

# ------
# Look at surveys data 
library(plyr)

surveys <- ddply(nnpp4, .(BatGroupName, 
                        SiteID, Easting, Northing, Precision, SurveyID,
                        SurveyDateStart, SurveyDateEnd, surveytime_start_unix, surveytime_end_unix, TrapID, 
                        TrapEasting, TrapNorthing, week_of_year, month, year), 
                 summarise,
                 no_bats = max(BatNumber, na.rm = T),
                no_nats = length(SpeciesID[SpeciesID == 23]),
                no_pygs = length(SpeciesID[SpeciesID == 25]),
                no_nocs = length(SpeciesID[SpeciesID == 19]))
                 #)


# do some more visualisations 
# change the levels of month so that they're in chronological order 
surveys$month = factor(surveys$month, levels = c('Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct'))

ggplot(surveys, aes(year, col = no_bats)) +geom_histogram(bins = 10)
ggplot(surveys, aes(month, no_bats)) + geom_boxplot()
ggplot(surveys, aes(month)) + geom_histogram(stat = 'count')
ggplot(surveys, aes(Easting, Northing, size=no_nats)) + geom_point(aes(col = month)) +facet_wrap(~year)
ggplot(surveys, aes(Easting, Northing, size=no_pygs, col = month)) + geom_point() +facet_wrap(~year)
ggplot(surveys, aes(Easting, Northing, size=no_nocs, col = month)) + geom_point() +facet_wrap(~year)
# how many surveys each year?
table(surveys$year)


# ------------------
# subset to only P. nathusius
nat <- subset(nnpp4, LatinName == 'Pipistrellus nathusii')
head(nat)
hist(nat$week_of_year)
hist(nat$month)

#How many surveys caught nathusius?
length(unique(nat$SurveyID))

#visualise
ggplot(nat, aes(Easting, Northing, size=BatNumber, col = Age)) + geom_point() +facet_wrap(~year)

length(unique(nat$RingNumber))
# ----------------
# subset to only P. pygmaeus
pyg <- subset(nnpp, LatinName == 'Pipistrellus pygmaeus')
head(pyg)
hist(pyg$week_of_year)


#--- ringed bats only

rngd <- nnpp[!is.na(nnpp$RingNumber), ]
head(rngd)

#list of unique ring numbers
ring_no <- unique(rngd$RingNumber)

ggplot(rngd, aes(Easting, Northing, col = LatinName)) + geom_point() +facet_wrap(~year)
















