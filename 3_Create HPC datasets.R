library(plyr); library(raster); library(ggplot2); library(dplyr); library(R2jags); library(rjags); library(mgcv); library(rgdal)


setwd("C:/Users/b284718/Google Drive/MSc Project ideas/Nats Pips/Bat Conservation Trust")
setwd("~/Dropbox/PhD/BCT_internship/analyses/data")

# ------------------------------------ #
# This script creates the datasets which are used in the JAGS models

## If the wind data has already been added to the surveys/ traps dfs then go to the point after that in this script.
# ------------------------------------ #

nnpp <- read.csv('nnpp_2014_18.csv', na.strings = c("", "NA")) # fill in blanks with NAs
head(nnpp)
# replace NA values in batnumber column with 0
nnpp$BatNumber[is.na(nnpp$BatNumber)] <- 0

# make a column for the months where the earliest month in the df is 1 - ie April becomes month 1 
nnpp$month_scaled <- nnpp$month - 3

# create column for survey start time since sunset 
nnpp$survey_start_since_sunset <- nnpp$surveytime_start_unix - nnpp$sunset_unix

# renumber sites so that theyre consecytive numbers
sites <- sort(unique(nnpp$SiteID))
sites_renumbered <- 1:length(sites)

sites_ref <- data.frame(sites, sites_renumbered)

for (i in 1:nrow(nnpp)){
  
  site_i <- nnpp$SiteID[i]
  
  nnpp$SiteID_x[i] <- sites_ref$sites_renumbered[sites_ref$sites == site_i]
  
}


# -- make a dataframe with sites and the number of surveys done at each one 
sites <- ddply(nnpp, .(SiteID_x, SiteGridReference, Easting_1km, Northing_1km, BatGroupName), summarise,
               no_visits = length(unique(SurveyID)))
          

# -- make a dataframe for individual trap events
# this df is for the observation model, so includes weather vars thought to influence detection of bats

traps <- ddply(nnpp, .(BatGroupName, 
                       SiteID_x, Easting_1km, Northing_1km, Precision, SurveyID,
                       SurveyDateStart, SurveyDateEnd, surveytime_start_unix, surveytime_end_unix, 
                       surveylength, survey_start_since_sunset,
                       TrapID, TrapEasting, TrapNorthing, 
                       Lure, month_scaled, year, 
                       lcm_class_names, StartTemp_2),#, wind_speed, wind_dir_factor), 
               summarise,
               no_bats = max(BatNumber),
               no_nats = length(SpeciesID[!is.na(SpeciesID) & SpeciesID == 23]),
               no_pygs = length(SpeciesID[!is.na(SpeciesID) & SpeciesID == 25]),
               no_nocs = length(SpeciesID[!is.na(SpeciesID) & SpeciesID == 19]))



# --  make a dataframe of surveys 
# this is for the survey level part of the model so doesn't need survey length, temp and wind data as these are just for the observation model

surveys <- ddply(traps,.(BatGroupName, 
                         SiteID_x, Easting_1km, Northing_1km, longitude, latitude,
                         SurveyID, SurveyDateStart, SurveyDateEnd, 
                         surveytime_start_unix, surveytime_end_unix, 
                         week_of_year, month_scaled, year, 
                         lcm_class_names), 
                 summarise,
                 no_bats = sum(no_bats),
                 no_nats = sum(no_nats),
                 no_pygs = sum(no_pygs),
                 no_nocs = sum(no_nocs))

# ------------------------ #
# Add wind data 
# ------------------------ #
# Here the surveys df is used to extract the relevent wind speed/direction values from the big wind data csv
# The survey df is used as this has fewer rows and the loop written below takes a long time (there is likelya faster way to do the extraction, but I haven't worked it out)
## Wind csv was created in get_wind_data.R script
wind_data <- read.csv('all_wind_data.csv')

# add a column to surveys with the rounded up long/lat to 0.5, which is the resolution of the wind data
surveys$longitude_2 <- round_any(surveys$longitude, 0.5)
surveys$latitude_2 <- round_any(surveys$latitude, 0.5)

# make sure start date column is in date format
surveys$SurveyDateStart_2<- format(strptime(as.character(surveys$SurveyDateStart), format="%d-%b-%y"), "%Y-%m-%d")

# Run big loop to get the data for each survey - TAKES HOURS
surveys$wind_speed <- NA
surveys$wind_direction <- NA
for(i in 1:nrow(surveys)){
  if (i %% 20 == 0){
    cat(paste0("row: ", i, "\n")) # this prints every 50th iteration of the loop so that you know it's running
  }
  
  for (j in 1:nrow(wind_data)){
 
    surveys$wind_speed[i][surveys$longitude_2[i] == wind_data$lon[j] & 
                            surveys$latitude_2[i] == wind_data$lat[j] &
                            surveys$SurveyDateStart_2[i] == wind_data$time[j]] = wind_data$speed[j]
    
    surveys$wind_direction[i][surveys$longitude_2[i] == wind_data$lon[j] & 
                                surveys$latitude_2[i] == wind_data$lat[j] &
                                surveys$SurveyDateStart_2[i] == wind_data$time[j]] = wind_data$dir[j]
  }
}

# make wind direction a factor
# N = 337.5:360 and 0:22.5
# NE = 22.5:67.5
# E = 67.5:112.5
# SE = 112,5:157.5
# S = 157.5:202.5
# SW = 202.5:247.5
# W =247.5:292.5
# NW = 292.5:337.5
library(dplyr)
surveys$wind_dir_factor <- case_when(surveys$wind_direction >22.4 & surveys$wind_direction <67.5 ~ 'NE',
                                     surveys$wind_direction > 67.4 & surveys$wind_direction < 112.5 ~ 'E',
                                     surveys$wind_direction >112.4 & surveys$wind_direction < 157.5 ~ 'SE', 
                                     surveys$wind_direction > 157.4 & surveys$wind_direction < 202.5 ~ 'S',
                                     surveys$wind_direction > 202.4 & surveys$wind_direction < 247.5 ~ 'SW',
                                     surveys$wind_direction > 247.4 & surveys$wind_direction < 292.5 ~ 'W', 
                                     surveys$wind_direction > 292.4  & surveys$wind_direction < 337.5 ~ 'NW',
                                     surveys$wind_direction > 337.4 & surveys$wind_direction < 360 ~ 'N',
                                     surveys$wind_direction > 0 & surveys$wind_direction < 22.5 ~ 'N')


# to add to traps data
surv_wind <- cbind(surveys$SurveyID, surveys$wind_dir_factor, surveys$wind_direction, surveys$wind_speed)
colnames(surv_wind) <- c('SurveyID', 'wind_dir_factor', 'wind_direction','wind_speed' )

traps_2 <- merge(traps, surv_wind, by = 'SurveyID')
write.csv(traps_2, 'traps.csv')
write.csv(surveys, 'surveys.csv')
# -----------------------------------------------------------#
# Start here if above steps have been completed previously #
# -----------------------------------------------------------#

## Read in the traps and surveys df if starting here
traps <- read.csv('traps.csv')
surveys <- read.csv('surveys.csv')

# make a column which tells you how often the site has been surveyed previously
site_id <- unique(surveys$SiteID_x)
surveys2 <- data.frame() 

for (i in 1:length(site_id)){
  
  site_i <- surveys[surveys$SiteID_x == site_id[i], ] # gets all the rows for each site 
  
  for (j in 1:nrow(site_i)){
    
    site_i$survey_event_id[j] = j # basically counting each row in the site specific df
    
  }
  surveys2 <- rbind(surveys2, site_i)
}

# now in the traps df, make a column which tells you the number each trap is within the survey (e.g 1, 2 or 3)
# also add in the survey event id 
survey_id <- unique(traps$SurveyID)
traps2 <- data.frame()

for (i in 1:length(survey_id)){
  
  survey_i <- traps[traps$SurveyID == survey_id[i], ] # gets all the rows for a specific survey
  survey_i$survey_event_id = surveys2$survey_event_id[surveys2$SurveyID == survey_i$SurveyID]
  
  for (j in 1:nrow(survey_i)){
    
    survey_i$trap_event_id[j] = j # basically counts the row number 
    
  }
  traps2 <- rbind(traps2, survey_i)
}
head(traps2)

# do some double checking of the data frames just created by using a random site  
# look at the counts of number of bats counted per survey and per trap 
sites[sites$SiteID_x == 121, ] # look at how many visits
surveys2[surveys2$SiteID_x == 121, ] # are there the same number of rows as visits in the site df
traps2[traps2$SiteID_x == 121, ] # are the trapIDs the same as in the above dataframe

# ---
# In order for the models to fit well, we need to remove some sites which are far away from the majority
# --- 

gb_outline <- readOGR(dsn = ".", layer = "GB_coarse_osgb")
plot(gb_outline)
points(sites$Easting_1km,sites$Northing_1km)
abline(h=6.5e+05)
abline(v=3e+05)

ggplot(sites, aes(x=Easting_1km, y = Northing_1km, col = BatGroupName)) + geom_point()
# groups to remove - Avon, Dorset, Cornwall, Jersey and North East Scotland 
rm_groups <- sites[sites$BatGroupName == 'Avon Bat Group' | sites$BatGroupName == 'Dorset Bat Group'| sites$BatGroupName == 'Cornwall Bat Group'|
                   sites$BatGroupName == 'Jersey Bat Group'| sites$BatGroupName == 'North East Scotland Bat Group',]

traps3 <- traps2[!traps2$BatGroupName %in% rm_groups$BatGroupName, ]
traps3 <- traps3[traps3$Easting_1km > 0, ] # 3 surveys with easting/northing as zero - have trap lcations though so will ammend these at a later date
ggplot(traps3, aes(x=Easting_1km, y = Northing_1km, col = BatGroupName)) + geom_point()
surveys3 <- surveys2[surveys2$SurveyID %in% traps3$SurveyID, ]

# standardise continuous variables 
## for the observation mmodel data
traps3$Easting_1km <- (traps3$Easting_1km-mean(traps3$Easting_1km))/sd(traps3$Easting_1km)
traps3$Northing_1km <- (traps3$Northing_1km-mean(traps3$Northing_1km))/sd(traps3$Northing_1km)
traps3$StartTemp_2 <- (traps3$StartTemp_2 - mean(na.omit(traps3$StartTemp_2)))/sd(na.omit(traps3$StartTemp_2))
traps3$surveylength <- (traps3$surveylength - mean(traps3$surveylength))/ sd(traps3$surveylength)
traps3$wind_speed <- (traps3$wind_speed - mean(na.omit(traps3$wind_speed)))/sd(na.omit(traps3$wind_speed))
traps3$survey_start_since_sunset <- (traps3$survey_start_since_sunset - mean(na.omit(traps3$survey_start_since_sunset)))/sd(na.omit(traps3$survey_start_since_sunset))
  
## for survey model 
surveys3$Easting_1km <- (surveys3$Easting_1km-mean(surveys3$Easting_1km))/sd(surveys3$Easting_1km)
surveys3$Northing_1km <- (surveys3$Northing_1km-mean(surveys3$Northing_1km))/sd(surveys3$Northing_1km)

# Check what the levels are for the categorical variables
# When these are put into the matrices below, factors are changed into their internal numeric values 
levels(traps3$wind_dir_factor)
levels(traps3$Lure)
levels(surveys3$lcm_class_names)

# make some matrices of trap counts and survey counts
#
obs_data <- data.frame(matrix(NA, nrow = nrow(traps3), ncol = 16))
names(obs_data) <- c('count', 'east', 'north',  
                     'month', 'year', 
                     'site_id', 'survey_id','survey_event_id', 
                     'trap_id',  'trap_event_id', 
                     'survey_dur', 'start_temp',
                     'wind_speed', 'wind_direction', 
                     'start_time_since_sunset', 'lure')
obs_data[, 1:16] <- c(traps3$no_nats, traps3$Easting_1km, traps3$Northing_1km, 
                      traps3$month_scaled, traps3$year,
                      traps3$SiteID_x, traps3$SurveyID, traps3$survey_event_id,
                      traps3$TrapID, traps3$trap_event_id, 
                      traps3$surveylength, traps3$StartTemp_2,#,
                       traps3$wind_speed, traps3$wind_dir_factor,
                      traps3$survey_start_since_sunset, traps3$Lure)
head(obs_data)

site_obsdata <- data.frame(matrix(NA,nrow=nrow(surveys3),ncol=9))
names(site_obsdata) <- c("count","east","north",
                         "month", "year",
                         "site_id", "survey_id", "survey_event_id", 
                         'lcm_class')
site_obsdata[, 1:9] <- c(surveys3$no_nats, surveys3$Easting_1km, surveys3$Northing_1km, 
                          surveys3$month_scaled, surveys3$year,
                          surveys3$SiteID_x, surveys3$SurveyID, surveys3$survey_event_id, 
                          surveys3$lcm_class_names)

head(site_obsdata)                         

plot(site_obsdata$east,site_obsdata$north)

# order by site_id
site_obsdata <- site_obsdata[order(site_obsdata$site_id),]
obs_data <- obs_data[order(obs_data$site_id),]

# rescale year in site obs data - this will be fitted as a random effect 
site_obsdata$year_rescaled <- (site_obsdata$year - 2013)

write.csv(site_obsdata,"site_obsdata_clipped.csv")
write.csv(obs_data,"obs_data_clipped.csv")

