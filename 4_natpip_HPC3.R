 library(raster);library(R2jags); library(rjags); library(mgcv); library(rgdal)

p <-  as.numeric(Sys.getenv("LSB_JOBINDEX")) # this gets the job number

#These are the options for k1 and k2 that we'll run:
k1s <- c(seq(5,15,1),seq(5,15,1),seq(5,15,1))
k2s <- c(rep(3,length(k1s)/3),rep(4,length(k1s)/3),rep(5,length(k1s)/3))


k1 <- k1s[p] #number of knots on each dimension of spatial grid
k2 <- k2s[p]
 
setwd("/gpfs/scratch/b284718/natpips")
#load data obkects
site_obsdata <- read.csv("site_obsdata_clipped.csv",row.names="X")
obs_data <- read.csv("obs_data_clipped.csv",row.names="X")
head(site_obsdata)
head(obs_data) 

# don't have wind data for pre 2015 from rWind package 
#remove surveys in 2014
site_obsdata <- site_obsdata[site_obsdata$year > 2014, ]
obs_data <- obs_data[obs_data$year > 2014, ]

# edit rescaled yrs
site_obsdata$year_rescaled <- site_obsdata$year_rescaled -1

#combine first and second months due to sparse data at start
#site_obsdata$month <- site_obsdata$month -1 
#obs_data$month <- obs_data$month -1 
#site_obsdata$month[site_obsdata$month==0]<-1 
#obs_data$month[obs_data$month==0] <-1 
# get the unuque years - in this case using the rescaled years
years <- site_obsdata$year_rescaled
maxyear <- max(years)

#The final object we must send is the 'key' vector created above which matches between rows of the trap and site level datasts
# This allows the JAGS code to identify which bit of the spline contributes to which trap-level observation
key <- vapply(1:nrow(obs_data),FUN=function(x) which(site_obsdata$site_id==obs_data$site_id[x] & site_obsdata$survey_event_id == obs_data$survey_event_id[x]),
              FUN.VALUE=1)

#Above matches each row in the big obs_data object to rows in t he site_obsdata object

# Now set up a background grid for interpolation
#gb_outline <- readOGR(dsn = ".", layer = "GB_coarse_osgb")
zone <- raster(matrix(0,nrow=50,ncol=50))
extent(zone) <- c(range(site_obsdata$east),range(site_obsdata$north))
plot(zone)
xy <- data.frame(xyFromCell(zone,1:ncell(zone)))

newdata1 <- data.frame(matrix(NA,nrow=nrow(xy),ncol=ncol(site_obsdata)))
names(newdata1) <- names(site_obsdata)
newdata1$count <- rep(0,nrow(xy))
newdata1$east <- xy$x
newdata1$north <- xy$y
newdata1$month<-1
newdata <- newdata1
months <- max(site_obsdata$month)

for(k in 1:months){
	newdata1$month<-k
	newdata <- rbind(newdata,newdata1)
}


## ------------------------------------------------------------------##

# -- Now fit the model --#

# First prepare the spline basis using jagam:

#Now repare splines using site-level data matrix (not trap level!)
# This is because the splines represent the big spatiotemporal pattern - variation between sites and months, not between traps
jagdata <- rbind(site_obsdata,newdata) # add the background grid points to the observed data

jags.ready <- jagam(count ~ te(east, north, bs="tp", k = k1) +
                      te(east, north, month, bs="tp", k = k2),	
                    data = jagdata, 
                    family = "poisson",
                    file ="jagam_simpleF.bug")

# Extract necessary objects from the jagam output:
S1<-jags.ready$jags.data$S1
S2<-jags.ready$jags.data$S2
X <- jags.ready$jags.data$X
zero <- jags.ready$jags.data$zero

#Now prepare the data objects that will be sent to the JAGS model:

kj1<- k1^2 # square of k1 - used in JAGS code for priors matrix dimensions
kj2<- length(zero)-kj1 # #also used for priors matrix dimensions
n1 <- nrow(site_obsdata) #number of site/month replicates for splines
n2 <- nrow(obs_data) #number of trap nigjht replicates for observation model

count <- obs_data$count #these are the observed counts
wind_speed <- obs_data$wind_speed
wind_direction <- obs_data$wind_direction
start_temp <- obs_data$start_temp
survey_dur <- obs_data$survey_dur
start_time_since_sunset <- obs_data$start_time_since_sunset
lure <- obs_data$lure
lcm_class <- site_obsdata$lcm_class

max_lcm_class <- max(lcm_class)
lcm_class.fac <- unique(lcm_class)
max_wind_direction <- max(na.omit(wind_direction))
wind_direction.fac <- unique(wind_direction)
max_lure <- max(lure)
lure_fac <- unique(lure)

# Now write the model file:
sink("Batgam_SPATEMP.txt")			 		
cat("
    model {
    #hyper-priors for mean presence spline
    for (i in 1:7) {
    lambda[i] ~ dgamma(.05,.005)
    }
    K1 <- S1[1:(kj1-1),1:(kj1-1)] * lambda[1]  + S1[1:(kj1-1),kj1:(2*(kj1-1))] * lambda[2] + S1[1:(kj1-1),(1+2*(kj1-1)):(3*(kj1-1))] * lambda[3] 
    K2 <- S2[1:(kj2),1:(kj2)] * lambda[4]  + S2[1:(kj2),(kj2+1):(2*(kj2))] * lambda[5] + S2[1:(kj2),(2*kj2+1):(3*(kj2))] * lambda[6] + S2[1:(kj2),(3*kj2+1):(4*(kj2))] * lambda[7] 
    
    # Priors on mean presence spline
    for(i in 1:1){b[i] ~ dnorm(0,0.0092)} 
    b[2:kj1] ~ dmnorm(zero[2:kj1],K1) 
    b[(kj1+1):(kj1+kj2)] ~ dmnorm(zero[(kj1+1):(kj1+kj2)],K2) 
    
    eta ~ dnorm(0,0.001) #random intercept on detection prob
    
    mu.beta[1:n1] <- X[1:n1,] %*% b ## linear predictor

 # priors on year random effects
  for(i in 1:maxyear){
  year.rand[i] ~ dnorm(0,0.0001)
  }
  
  # priors on categorical variables
  lcm_class.fac[1] <- 0 
  for(i in 2:max_lcm_class) {
    lcm_class.fac[i] ~ dnorm(0,0.0001) 
  }
    
  wind_direction.fac[1] <- 0 
  for(i in 2:max_wind_direction){
  wind_direction.fac[i] ~ dnorm(0, 0.0001)
  }

lure_fac[1] <- 0
for(i in 1:max(lure)){
lure_fac[i] ~ dnorm(0, 0.0001)
}

 #priors for betas 1-4
 beta1 ~ dnorm(0,0.0001)
 beta2 ~ dnorm(0,0.0001)
 beta3 ~ dnorm(0,0.0001)
 beta4 ~ dnorm(0,0.0001)

    for (i in 1:n1) { 
    log(mu[i]) <-  mu.beta[i] + year.rand[years[i]] + lcm_class.fac[lcm_class[i]]
    z[i] ~ dpois(mu[i]) # latent true presence vector
    }

    for (i in 1:n2) { 
    count[i] ~ dbin(p[i],z[key[i]]) ## observation model - z is true abund p is detection prob
    start_temp[i] ~ dnorm(0,0.001) ## priors on predictors with missing values
    wind_speed[i] ~ dnorm(0, 0.001)
    logit(p[i]) <- eta + beta1 *  wind_speed[i] + wind_direction.fac[wind_direction[i]] + beta2 * start_temp[i] + beta3 * survey_dur[i] + beta4 * start_time_since_sunset[i] + lure_fac[lure[i]] # fixed value of detection prob on logit scale
    }

    }
    ", fill=TRUE)
sink()
 
# Initial values
inits <- function(){ list(
  z=site_obsdata$count, #including this helps initial model fitting
	b=jags.ready$jags.ini$b,
	year.rand=rep(0,maxyear),
  
  lambda=jags.ready$jags.ini$lambda,
  eta=rnorm(1,0,0.1)
  #
)}


n.iter <- 10000
n.burn <-5000
n.thin <- 10
n.chains<-3
params <-c("b","eta","lambda","beta1","beta2", "beta3", "beta4",
           "wind_direction.fac", "lcm_class.fac", "year.rand", "lure_fac")#,


# Fit the model - might take a while!
model.fit <- do.call(jags, list(c("count", "n1", "n2","X", "S1", "S2","zero", "kj1", "kj2", "key", 
                                  "years", "maxyear", "lcm_class", "max_lcm_class",
                                  "wind_speed", "wind_direction", "max_wind_direction", "start_temp", 
                                  "survey_dur", "start_time_since_sunset", "lure", "max_lure"),
                     inits, params,"Batgam_SPATEMP.txt",    
                                n.chains, n.iter, n.burn, n.thin))

save.image(file=paste("Natpip_clip2_k1_",k1,"_k2_",k2,"beta20005",".RData",sep=""))
	