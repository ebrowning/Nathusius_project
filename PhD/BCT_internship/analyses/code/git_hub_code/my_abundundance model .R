library(plyr); library(raster); library(ggplot2); library(dplyr); library(R2jags); library(rjags); library(mgcv); library(coda)


setwd("~/Dropbox/PhD/BCT_internship/analyses/data/")

nnpp <- read.csv('nnpp_2014_18.csv', na.strings = c("", "NA")) # fill in blanks with NAs
head(nnpp)
# replace NA values in batnumber column with 0
nnpp$BatNumber[is.na(nnpp$BatNumber)] <- 0

# make a column for the months where the earliest month in the df is 1 - ie March becomes month 1 
nnpp$month_scaled <- nnpp$month - 2

# make a column for season based on week of the year (as migrants likely to arrive back in the UK in mid august, makes more sense than roughly splitting by month)
# 1 = spring [week 11 - week 21]; 2 = summer [week 22 - week 32] ; 3 = autumn [week 33 - 44]

nnpp$season = case_when(nnpp$week_of_year >10 & nnpp$week_of_year <22 ~ 1,
                        nnpp$week_of_year > 21 & nnpp$week_of_year < 32 ~ 2,
                        nnpp$week_of_year >33 ~ 3)


# renumber sites so that theyre consecytive numbers
sites <- sort(unique(nnpp$SiteID))
sites_renumbered <- 1:length(sites)

sites_ref <- data.frame(sites, sites_renumbered)

for (i in 1:nrow(nnpp)){
  
  site_i <- nnpp$SiteID[i]
  
  nnpp$SiteID_x[i] <- sites_ref$sites_renumbered[sites_ref$sites == site_i]
  
}


# -- make a dataframe with sites and the number of surveys done at each one 
sites <- ddply(nnpp, .(SiteID_x, SiteGridReference, Easting_1km, Northing_1km), summarise,
               no_visits = length(unique(SurveyID)))
          

# -- make a dataframe for individual trap events

traps <- ddply(nnpp, .(BatGroupName, 
                       SiteID_x, Easting_1km, Northing_1km, Precision, SurveyID,
                       SurveyDateStart, SurveyDateEnd, surveytime_start_unix, surveytime_end_unix, 
                       TrapID, TrapEasting, TrapNorthing, 
                       week_of_year, month_scaled, season, year), 
               summarise,
               no_bats = max(BatNumber),
               no_nats = length(SpeciesID[SpeciesID == 23]),
               no_pygs = length(SpeciesID[SpeciesID == 25]),
               no_nocs = length(SpeciesID[SpeciesID == 19]))



# --  make a dataframe of surveys 

surveys <- ddply(traps,.(BatGroupName, 
                         SiteID_x, Easting_1km, Northing_1km, Precision, SurveyID,
                         SurveyDateStart, SurveyDateEnd, surveytime_start_unix, surveytime_end_unix, 
                         week_of_year, month_scaled, season, year), 
                 summarise,
                 no_bats = sum(no_bats),
                 no_nats = sum(no_nats),
                 no_pygs = sum(no_pygs),
                 no_nocs = sum(no_nocs))

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


# do some double checking of the data frames just created by using a random site  
# look at the counts of number of bats counted per survey and per trap 
sites[sites$SiteID_x == 201, ] # look at how many visits
surveys2[surveys2$SiteID_x == 201, ] # are there the same number of rows as visits in the site df
traps2[traps2$SiteID_x == 201, ] # are the trapIDs the same as in the above dataframe
nnpp[nnpp$SiteID_x == 201, ] # check back with original dataframe to make sure the counts of bats are right

# make some matrices of trap counts and survey counts
obs_data <- data.frame(matrix(NA, nrow = nrow(traps), ncol = 9))
names(obs_data) <- c('nats_count', 'site_id', 'survey_id', 'survey_event_id', 'trap_id', 'trap_event_id', 'month', 'east', 'north')
obs_data[, 1:9] <- c(traps2$no_nats, traps2$SiteID_x, traps2$SurveyID, traps2$survey_event_id, traps2$TrapID, traps2$trap_event_id, traps2$month_scaled, traps2$Easting_1km, traps2$Northing_1km)
head(obs_data)

site_obsdata <- data.frame(matrix(NA,nrow=nrow(surveys2),ncol=7))
names(site_obsdata) <- c("count","east","north","month","site_id", "survey_id", "survey_event_id")
site_obsdata[, 1:7] <- c(surveys2$no_nats, surveys2$Easting_1km, surveys2$Northing_1km, surveys2$month_scaled, surveys2$SiteID_x, surveys2$SurveyID, 
                         surveys2$survey_event_id)
head(site_obsdata)                         

# make a vector of each unique season 
month = unique(nnpp$month_scaled)

month#The final object we must send is the 'key' vector created above which matches between rows of the trap and site level datasts
# This allows the JAGS code to identify which bit of the spline contributes to which trap-level observation
key <- vapply(1:nrow(obs_data),FUN=function(x) which(site_obsdata$site_id==obs_data$site_id[x] & site_obsdata$survey_event_id == obs_data$survey_event_id[x]),
              FUN.VALUE=1)
#Above matches each row in the big obs_data object to rows in t he site_obsdata object

## ------------------------------------------------------------------##

# -- Now fit the model --#

# First prepare the spline basis using jagam:

#Number of knots for the gam splines:
k1 <- 9 #number of knots on each dimension of spatial grid
k2 <-  7#seasons-1 #number of knots for temporal interaction spline

#Now repare splines using site-level data matrix (not trap level!)
# This is because the splines represent the big spatiotemporal pattern - variation between sites and months, not between traps

jags.ready <- jagam(count ~ te(east, north, bs="tp", k = k1) +
                      te(east, north, month, bs="tp", k = k2),	
                    data = site_obsdata, 
                    family = "poisson",
                    file ="jagam_simpleF.bug")

# Extract necessary objects from the jagam output:
S1<-jags.ready$jags.data$S1
S2<-jags.ready$jags.data$S2
X <- jags.ready$jags.data$X
zero <- jags.ready$jags.data$zero
rm(jags.ready)

#Now prepare the data objects that will be sent to the JAGS model:

kj1<- k1^2 # square of k1 - used in JAGS code for priors matrix dimensions
kj2<- length(zero)-kj1 # #also used for priors matrix dimensions
n1 <- nrow(site_obsdata) #number of site/month replicates for splines
n2 <- nrow(obs_data) #number of trap nigjht replicates for observation model

count <- obs_data$nats_count #these are the observed counts

# Now write the model file:
sink("Batgam_SPATEMP.txt")			 		
cat("
    model {
    #hyper-priors for mean presence spline

	  for (i in 1:7) {
    lambda[i] ~ dgamma(0.05, 0.005)
    }
 		K1 <- S1[1:(kj1-1),1:(kj1-1)] * lambda[1]  + S1[1:(kj1-1),kj1:(2*(kj1-1))] * lambda[2] + S1[1:(kj1-1),(1+2*(kj1-1)):(3*(kj1-1))] * lambda[3] 
  	K2 <- S2[1:(kj2),1:(kj2)] * lambda[4]  + S2[1:(kj2),(kj2+1):(2*(kj2))] * lambda[5] + S2[1:(kj2),(2*kj2+1):(3*(kj2))] * lambda[6] + S2[1:(kj2),(3*kj2+1):(4*(kj2))] * lambda[7] 
 
 	# Priors on mean presence spline
  	for(i in 1:1){
  	  b[i] ~ dnorm(0, 0.0001)
  	  } 
   	b[2:kj1] ~ dmnorm(zero[2:kj1],K1) 
  	b[(kj1+1):(kj1+kj2)] ~ dmnorm(zero[(kj1+1):(kj1+kj2)],K2) 

		eta ~dnorm(0, 0.001) #random intercept on detection prob
  	
		mu.beta[1:n1] <- X[1:n1,] %*% b ## linear predictor

    for (i in 1:n1) { 
	  log(mu[i]) =  mu.beta[i]
    z[i] ~ dpois(mu[i]) # latent true presence vector
    }
    
 		for (i in 1:n2) { 
		count[i] ~ dbin(p[i],z[key[i]]) ## observation model - z is true abund p is detection prob
    logit(p[i]) <- eta # fixed value of detection prob on logit scale
}
    }
    ",fill=TRUE)
sink()
 
# Initial values
inits <- function(){ list(
  z=site_obsdata$count+10, #including this helps initial model fitting
  lambda=rep(0.005,7),
  eta=rnorm(1,0,0.1)
  #
)}


n.iter <- 10000
n.burn <-8000
n.thin <- 8
n.chains<-4
params <-c("b","eta","lambda")#,

# Fit the model - might take a while!
model.fit <- do.call(jags, list(c("count","n1","n2","X","S1","S2","zero","kj1","kj2","key"), 
                                inits, params,"Batgam_SPATEMP.txt",    
                                n.chains, n.iter, n.burn, n.thin))

# Check for convergence - ideally Rhat for all nodes <1.05 
hist(model.fit$BUGSoutput$summary[,"Rhat"])
sum(model.fit$BUGSoutput$summary[,"Rhat"]>1.05)

# Keep updating the model if necessary to get to convergence-
model.fit <- update(model.fit,n.iter=5000,n.thin=n.thin)

 # Explore model outputs:
# Posterior distribution of p_d 
p_d # recall true value of detect prob
p_d_samples <- 1/(1+exp(-(model.fit$BUGSoutput$sims.list$eta))) #extract posterior samples of eta, inverse logit transformed to p_d

hist(p_d_samples,breaks=30)

model.fit$BUGSoutput$DIC
model.fit$BUGSoutput$summary

model.fit.mcmc <- as.mcmc(model.fit)
summary(model.fit.mcmc)

densplot(model.fit.mcmc)


# Explore model-predicted patterns:  

b.samps <- model.fit$BUGSoutput$mean$b #get posterior means of the b matrix (ie the fitted spline basis)

est.counts <- exp(X%*%b.samps) # generate the model estimated true counts - matrix product of spline basis and X matrix generated by jagam, inverse log transformed 

# Now create out object to plot the estimated counts across space and time
maps <- data.frame(x=site_obsdata$east,y=site_obsdata$north, month = site_obsdata$month, est.counts)

colsc <- heat.colors(30)
bks <- seq(0,max(maps$est.counts),length.out=30) #color scheme for predicted counts
colsc[1] <- "black"
colsc2 <- rev(terrain.colors(50))
bks2 <- seq(0,max(cellStats(dists,max)),length.out=50) #color scheme for base distribution

# Make a plot for each month, compare model to base distribution
for(k in 1:months){
  plot(dists[[distpick[k]]], main=paste("Month",k),col=colsc2,breaks=bks2,legend=F)
  map1 <- subset(maps,month==k)
  points(map1$x,map1$y,pch=16,col=colsc[findInterval(map1$est.counts,bks)],cex=1.5)
  par(xpd=T)
  xleg <- seq(5,50,5)
  points(x= xleg,y=rep(-30,length(xleg)),pch=16,cex=2, col= colsc[floor(seq(1,length(colsc),length.out=length(xleg)))])
  text(x= c(5,50),y=rep(-30,2),c("0",round(max(bks),1)),pos=1)
  text(x= c(80),y=-30,c("Est. count"),pos=1)
}








